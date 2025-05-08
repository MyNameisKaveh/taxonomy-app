from flask import Flask, jsonify, request
import requests
import traceback
import wikipedia
from Bio import Entrez # اضافه کردن Entrez
import time

# --- کتابخانه‌های جدید ---
from googletrans import Translator, LANGUAGES
from langdetect import detect, LangDetectException

app = Flask(__name__)
GBIF_API_URL_MATCH = "https://api.gbif.org/v1/species/match"

# !!! ایمیل خودتون رو اینجا برای Entrez تنظیم کنید !!!
Entrez.email = "YOUR_ACTUAL_EMAIL@example.com" 

# --- تابع تشخیص زبان و ترجمه ---
def translate_to_english_if_needed(text):
    try:
        detected_lang = detect(text)
        print(f"[TRANSLATE] Detected language: {detected_lang} for '{text}'")
        if detected_lang == 'fa': # اگر فارسی بود
            translator = Translator()
            # ممکن است نیاز به تعیین صریح منبع و مقصد باشد
            translation = translator.translate(text, src='fa', dest='en')
            translated_text = translation.text
            print(f"[TRANSLATE] Translated '{text}' (fa) to '{translated_text}' (en)")
            return translated_text, 'fa' # نام ترجمه شده و زبان اصلی
        else:
            # اگر انگلیسی یا زبان دیگری بود، خود متن رو برگردون
            return text, detected_lang
    except LangDetectException:
        print(f"[TRANSLATE] Language detection failed for '{text}'. Assuming English.")
        return text, 'unknown' # اگر تشخیص زبان شکست خورد
    except Exception as e:
        print(f"[TRANSLATE] Error during translation: {e}")
        return text, 'error' # در صورت خطای ترجمه، متن اصلی رو برگردون

# --- تابع جستجوی NCBI (همان کد کامل قبلی) ---
def get_best_ncbi_suggestion_flexible(common_name, max_ids_to_check=5):
    # ... (کد کامل و صحیح تابع get_best_ncbi_suggestion_flexible از پیام قبلی اینجا قرار می‌گیرد) ...
    print(f"\n--- Processing '{common_name}' with NCBI Entrez (Flexible Ranks) ---")
    scientific_name_suggestion = None
    try:
        search_term = f"{common_name}[Common Name] OR {common_name}[Organism]"
        handle = Entrez.esearch(db="taxonomy", term=search_term, retmax=max_ids_to_check, sort="relevance")
        record = Entrez.read(handle)
        handle.close()
        id_list = record["IdList"]

        if not id_list:
            time.sleep(0.34)
            handle = Entrez.esearch(db="taxonomy", term=common_name, retmax=max_ids_to_check, sort="relevance")
            record = Entrez.read(handle)
            handle.close()
            id_list = record["IdList"]
            if not id_list:
                print(f"  No TaxIDs found for '{common_name}' even after retry.")
                return None
        
        print(f"  Found {len(id_list)} potential TaxIDs: {id_list}")
        ids_to_fetch = id_list[:max_ids_to_check] 
        if not ids_to_fetch: return None

        wait_time = 0.2 + (0.34 * len(ids_to_fetch)) 
        time.sleep(wait_time) 
        
        summary_handle = Entrez.esummary(db="taxonomy", id=",".join(ids_to_fetch), retmode="xml")
        summary_records = Entrez.read(summary_handle)
        summary_handle.close()
        
        candidates = {"species": None, "subspecies": None, "genus": None, "family": None, "other_valid": None}
        
        for summary_record in summary_records:
            sci_name = summary_record.get("ScientificName")
            rank_ncbi = summary_record.get("Rank", "N/A").lower()
            tax_id_processed = summary_record.get("Id") 
            print(f"    - Candidate TaxID {tax_id_processed}: '{sci_name}', Rank: {rank_ncbi}")
            if sci_name:
                if rank_ncbi == "species":
                    if not candidates["species"]: candidates["species"] = sci_name; print(f"      -> Found Species Match: '{sci_name}'")
                elif rank_ncbi == "subspecies":
                    if not candidates["subspecies"]: candidates["subspecies"] = sci_name; print(f"      -> Found Subspecies Match: '{sci_name}'")
                elif rank_ncbi == "genus":
                     if not candidates["genus"]: candidates["genus"] = sci_name; print(f"      -> Found Genus Match: '{sci_name}'")
                elif rank_ncbi == "family":
                     if not candidates["family"]: candidates["family"] = sci_name; print(f"      -> Found Family Match: '{sci_name}'")
                elif not candidates["other_valid"]: candidates["other_valid"] = sci_name; print(f"      -> Found Other Valid Scientific Name: '{sci_name}' (Rank: {rank_ncbi})")
        
        if candidates["species"]: scientific_name_suggestion = candidates["species"]
        elif candidates["subspecies"]:
            parts = candidates["subspecies"].split(); scientific_name_suggestion = " ".join(parts[:2]) if len(parts) >= 2 else candidates["subspecies"]
            print(f"  -> Using Species name from Subspecies: '{scientific_name_suggestion}'")
        elif candidates["genus"]: scientific_name_suggestion = candidates["genus"]; print(f"  -> Using Genus name as best match: '{candidates['genus']}'")
        elif candidates["family"]: scientific_name_suggestion = candidates["family"]; print(f"  -> Using Family name as best match: '{candidates['family']}'")
        elif candidates["other_valid"]: scientific_name_suggestion = candidates["other_valid"]; print(f"  -> Using first other valid name as best match: '{candidates['other_valid']}'")
        else: print(f"  No suitable scientific name found according to prioritization.")

    except Exception as e: print(f"  Error querying NCBI Entrez: {type(e).__name__} - {e}")
    return scientific_name_suggestion

# --- Endpoint جدید برای راهنما ---
@app.route('/api/suggest_name', methods=['GET'])
def suggest_name_endpoint():
    query = request.args.get('query')
    if not query:
        return jsonify({"error": "Query parameter is required"}), 400

    # تشخیص زبان و ترجمه در صورت نیاز
    name_to_search, original_lang = translate_to_english_if_needed(query)

    # جستجو با نام انگلیسی (اصلی یا ترجمه شده)
    suggestion = get_best_ncbi_suggestion_flexible(name_to_search)

    response_data = {
        "query": query,
        "searched_term_en": name_to_search,
        "original_language": original_lang,
        "scientific_name_suggestion": suggestion
    }

    if suggestion:
        return jsonify(response_data)
    else:
        response_data["message"] = f"No definitive scientific name suggestion found for '{query}'."
        return jsonify(response_data), 404

# --- تابع اصلی جستجوی طبقه‌بندی (main_handler) و تابع get_wikipedia_image_url ---
# این توابع باید مثل قبل در فایل باشند، اما دیگر نیازی به تغییر ندارند
# (کد کامل آنها در پیام‌های خیلی قبل‌تر هست)
@app.route('/', defaults={'path': ''}, methods=['GET', 'POST', 'OPTIONS'])
@app.route('/<path:path>', methods=['GET', 'POST', 'OPTIONS'])
def main_handler(path=None):
    # ... (کد کامل main_handler که تصویر را از ویکی‌پدیا می‌گیرد و به GBIF وصل می‌شود) ...
    # این تابع بدون تغییر باقی می‌ماند
    common_headers = {'Access-Control-Allow-Origin': '*'}
    if request.method == 'OPTIONS':
        cors_headers = {**common_headers, 'Access-Control-Allow-Methods': 'GET, POST, OPTIONS', 'Access-Control-Allow-Headers': 'Content-Type, Authorization', 'Access-Control-Max-Age': '3600'}
        return ('', 204, cors_headers)

    species_name_query_original = ""
    if request.method == 'GET': species_name_query_original = request.args.get('name')
    # ... (گرفتن از POST)

    if not species_name_query_original: return jsonify({"error": "Parameter 'name' required"}), 400

    # اینجا دیگر نیازی به دیکشنری تبدیل نام رایج نیست
    # چون جستجوی اصلی فقط با نام علمی دقیق انجام می‌شود
    species_name_for_gbif = species_name_query_original # فقط نامی که کاربر داده
    print(f"[MAIN_HANDLER] Processing request for: '{species_name_for_gbif}'")

    params_gbif = {"name": species_name_for_gbif, "verbose": "true"}
    classification_data = {"searchedName": species_name_query_original}
    gbif_error_message = None
    gbif_scientific_name = None

    # ... (بلاک try/except کامل برای ارتباط با GBIF مثل قبل) ...
    try:
        api_response_gbif = requests.get(GBIF_API_URL_MATCH, params=params_gbif, timeout=10)
        api_response_gbif.raise_for_status()
        data_gbif = api_response_gbif.json()
        # ... (منطق پردازش پاسخ GBIF) ...
        if not data_gbif or data_gbif.get("matchType") == "NONE" or data_gbif.get("confidence", 0) < 30:
            gbif_error_message = f"موجودی با نام '{species_name_for_gbif}' در GBIF پیدا نشد یا نتیجه با اطمینان کافی نبود."
            classification_data["matchType"] = data_gbif.get("matchType", "NONE"); classification_data["confidence"] = data_gbif.get("confidence")
        else:
            gbif_scientific_name = data_gbif.get("scientificName")
            classification_data.update({
                "scientificName": gbif_scientific_name, "kingdom": data_gbif.get("kingdom"), "phylum": data_gbif.get("phylum"),
                "class": data_gbif.get("class"), "order": data_gbif.get("order"), "family": data_gbif.get("family"),
                "genus": data_gbif.get("genus"), "species": data_gbif.get("species") if data_gbif.get("speciesKey") and data_gbif.get("species") else None,
                "usageKey": data_gbif.get("usageKey"), "confidence": data_gbif.get("confidence"), "matchType": data_gbif.get("matchType"),
                "status": data_gbif.get("status"), "rank": data_gbif.get("rank")
            })   
    except requests.exceptions.Timeout: gbif_error_message = "خطا: زمان پاسخگویی GBIF طول کشید."; print(f"[GBIF_ERR] Timeout for: {species_name_for_gbif}")
    except requests.exceptions.HTTPError as http_err_gbif: gbif_error_message = f"خطای سرور GBIF: {http_err_gbif}"; print(f"[GBIF_ERR] HTTPError: {gbif_error_message}")
    except requests.exceptions.RequestException as e_gbif_req: gbif_error_message = f"خطای شبکه GBIF: {str(e_gbif_req)}"; print(f"[GBIF_ERR] RequestException: {str(e_gbif_req)}")
    except Exception as e_gbif_generic: gbif_error_message = "خطای داخلی GBIF."; print(f"[GBIF_ERR] Generic Error: {str(e_gbif_generic)}"); traceback.print_exc()

    # تابع get_wikipedia_image_url بدون تغییر باقی می ماند
    # ... (کد کامل تابع get_wikipedia_image_url از پیام‌های قبلی اینجا قرار می‌گیرد) ...
    def get_wikipedia_image_url(species_name_from_user, scientific_name_from_gbif=None):
        # ... (همان کد کامل و دقیق شده قبلی) ...
        search_candidates=[]; clean_scientific_name=None; clean_scientific_name_for_filename=None
        if scientific_name_from_gbif: temp_clean_name=scientific_name_from_gbif.split('(')[0].strip();
        if temp_clean_name: clean_scientific_name=temp_clean_name; search_candidates.append(clean_scientific_name); clean_scientific_name_for_filename=clean_scientific_name.lower().replace(" ","_")
        user_name_for_filename=None
        if species_name_from_user: should_add_user_name=True
        if clean_scientific_name:
            if species_name_from_user.lower()==clean_scientific_name.lower(): should_add_user_name=False
        if should_add_user_name: search_candidates.append(species_name_from_user)
        user_name_for_filename=species_name_from_user.lower().replace(" ","_")
        if not search_candidates: print("[WIKI_IMG] No search terms provided."); return None
        print(f"[WIKI_IMG] Attempting image for candidates: {search_candidates}"); wikipedia.set_lang("en")
        avoid_keywords_in_filename=["map","range","distribution","locator","chart","diagram","logo","icon","disambig","sound","audio","timeline","scale","reconstruction","skeleton","skull","footprint","tracks","scat","phylogeny","cladogram","taxonomy","taxobox"]
        priority_keywords=[]
        if clean_scientific_name_for_filename: priority_keywords.append(clean_scientific_name_for_filename)
        if "_" in clean_scientific_name_for_filename: priority_keywords.append(clean_scientific_name_for_filename.split("_")[0])
        if user_name_for_filename:
            if not (clean_scientific_name_for_filename and user_name_for_filename==clean_scientific_name_for_filename): priority_keywords.append(user_name_for_filename)
            if "_" in user_name_for_filename: user_genus_equivalent=user_name_for_filename.split("_")[0]
            add_user_genus=True
            if clean_scientific_name_for_filename and "_" in clean_scientific_name_for_filename:
                if user_genus_equivalent==clean_scientific_name_for_filename.split("_")[0]: add_user_genus=False
            elif clean_scientific_name_for_filename and user_genus_equivalent==clean_scientific_name_for_filename : add_user_genus=False
            if add_user_genus: priority_keywords.append(user_genus_equivalent)
        priority_keywords=list(filter(None,dict.fromkeys(priority_keywords))); print(f"[WIKI_IMG] Priority keywords for image filename: {priority_keywords}")
        processed_search_terms=set()
        while search_candidates: term_to_search=search_candidates.pop(0)
        if not term_to_search or term_to_search in processed_search_terms: continue
        processed_search_terms.add(term_to_search); print(f"[WIKI_IMG] Trying Wikipedia search for term: '{term_to_search}'")
        try: wiki_page=None
        try: wiki_page=wikipedia.page(term_to_search,auto_suggest=True,redirect=True); print(f"[WIKI_IMG] Found page directly: '{wiki_page.title}' for '{term_to_search}'")
        except wikipedia.exceptions.PageError: print(f"[WIKI_IMG] Page not found directly for '{term_to_search}'. Trying wikipedia.search()."); search_results=wikipedia.search(term_to_search,results=1)
        if not search_results: print(f"[WIKI_IMG] No search results in Wikipedia for: '{term_to_search}'"); continue
        page_title_to_get=search_results[0]
        if page_title_to_get in processed_search_terms: print(f"[WIKI_IMG] Page '{page_title_to_get}' already processed, skipping."); continue
        print(f"[WIKI_IMG] Found page via search: '{page_title_to_get}' for '{term_to_search}'"); wiki_page=wikipedia.page(page_title_to_get,auto_suggest=False,redirect=True)
        if wiki_page and wiki_page.images: print(f"[WIKI_IMG] Page '{wiki_page.title}' images (up to 5): {wiki_page.images[:5]}"); candidate_images_with_scores=[]
        for img_url in wiki_page.images: img_url_lower=img_url.lower()
        if not any(ext in img_url_lower for ext in['.png','.jpg','.jpeg']): continue
        if any(keyword in img_url_lower for keyword in avoid_keywords_in_filename): continue
        score=0
        for pk_word in priority_keywords:
            if pk_word in img_url_lower: score+=5; filename_part=img_url_lower.split('/')[-1]
            if filename_part.startswith(pk_word): score+=3
            if clean_scientific_name_for_filename and pk_word==clean_scientific_name_for_filename and pk_word in filename_part : score+=5
        if img_url_lower.endswith('.svg'): score-=1; candidate_images_with_scores.append({'url':img_url,'score':score})
        if not candidate_images_with_scores: print(f"[WIKI_IMG] No images passed filter for '{wiki_page.title}'"); continue
        sorted_images=sorted(candidate_images_with_scores,key=lambda x:x['score'],reverse=True); print(f"[WIKI_IMG] Sorted suitable images (top 3 with scores): {[{'url':i['url'][-50:],'score':i['score']} for i in sorted_images[:3]]}")
        if sorted_images and sorted_images[0]['score']>0: best_image_url=sorted_images[0]['url']
        if best_image_url.startswith("//"): best_image_url=":"+best_image_url
        print(f"[WIKI_IMG] Best image found: {best_image_url} with score {sorted_images[0]['score']}"); return best_image_url
        else: print(f"[WIKI_IMG] No image found with positive score for '{wiki_page.title}'.")
        elif wiki_page: print(f"[WIKI_IMG] No images listed on Wikipedia page: '{wiki_page.title}'")
        except wikipedia.exceptions.DisambiguationError as e: print(f"[WIKI_IMG] Disambiguation for '{term_to_search}'. Options: {e.options[:3]}")
        if e.options: new_candidate=e.options[0]
        if new_candidate not in processed_search_terms and new_candidate not in search_candidates: search_candidates.append(new_candidate); print(f"[WIKI_IMG] Added disambiguation option '{new_candidate}' to search candidates.")
        continue
        except wikipedia.exceptions.PageError: print(f"[WIKI_IMG] Wikipedia PageError for term: '{term_to_search}'"); continue
        except Exception as e: print(f"[WIKI_IMG] Generic error for '{term_to_search}': {str(e)}"); traceback.print_exc(); continue
        print(f"[WIKI_IMG] No suitable Wikipedia image URL after all attempts for initial candidates."); return None

    # --- بخش نهایی main_handler ---
    wiki_image_url = get_wikipedia_image_url(species_name_query_original, gbif_scientific_name) 
    if wiki_image_url: classification_data["imageUrl"] = wiki_image_url
    final_data = {k: v for k, v in classification_data.items() if v is not None}
    if gbif_error_message and not final_data.get("scientificName"): 
        if final_data.get("imageUrl"): final_data["gbifLookupMessage"] = gbif_error_message; return jsonify(final_data), 200, common_headers
        return jsonify({"message": gbif_error_message, "searchedName": species_name_query_original}), 404, common_headers
    return jsonify(final_data), 200, common_headers
