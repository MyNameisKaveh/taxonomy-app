from flask import Flask, jsonify, request
import requests
import traceback
import wikipedia
from Bio import Entrez # اضافه شده
import time            # اضافه شده

app = Flask(__name__)
GBIF_API_URL_MATCH = "https://api.gbif.org/v1/species/match"

# !!! بسیار مهم: ایمیل واقعی خود را اینجا وارد کنید !!!
Entrez.email = "YOUR_ACTUAL_EMAIL@example.com"
# !!! بسیار مهم: ایمیل واقعی خود را اینجا وارد کنید !!!


# --- تابع جستجوی تصویر از ویکی‌پدیا (بدون تغییر) ---
def get_wikipedia_image_url(species_name_from_user, scientific_name_from_gbif=None):
    search_candidates = []
    clean_scientific_name = None 
    clean_scientific_name_for_filename = None

    if scientific_name_from_gbif:
        temp_clean_name = scientific_name_from_gbif.split('(')[0].strip()
        if temp_clean_name:
            clean_scientific_name = temp_clean_name
            search_candidates.append(clean_scientific_name)
            clean_scientific_name_for_filename = clean_scientific_name.lower().replace(" ", "_")
    
    user_name_for_filename = None
    if species_name_from_user:
        should_add_user_name = True
        if clean_scientific_name: 
            if species_name_from_user.lower() == clean_scientific_name.lower():
                should_add_user_name = False
        
        if should_add_user_name:
            search_candidates.append(species_name_from_user)
        
        user_name_for_filename = species_name_from_user.lower().replace(" ", "_")
    
    if not search_candidates:
        print("[WIKI_IMG] No search terms provided.")
        return None

    print(f"[WIKI_IMG] Attempting image for candidates: {search_candidates}")
    wikipedia.set_lang("en")

    avoid_keywords_in_filename = ["map", "range", "distribution", "locator", "chart", "diagram", "logo", "icon", "disambig", "sound", "audio", "timeline", "scale", "reconstruction", "skeleton", "skull", "footprint", "tracks", "scat", "phylogeny", "cladogram", "taxonomy", "taxobox"]
    
    priority_keywords = []
    if clean_scientific_name_for_filename: 
        priority_keywords.append(clean_scientific_name_for_filename)
        if "_" in clean_scientific_name_for_filename: 
            priority_keywords.append(clean_scientific_name_for_filename.split("_")[0])
            
    if user_name_for_filename:
        if not (clean_scientific_name_for_filename and user_name_for_filename == clean_scientific_name_for_filename):
            priority_keywords.append(user_name_for_filename)
        if "_" in user_name_for_filename:
            user_genus_equivalent = user_name_for_filename.split("_")[0]
            add_user_genus = True
            if clean_scientific_name_for_filename and "_" in clean_scientific_name_for_filename:
                 if user_genus_equivalent == clean_scientific_name_for_filename.split("_")[0]:
                     add_user_genus = False
            elif clean_scientific_name_for_filename and user_genus_equivalent == clean_scientific_name_for_filename : 
                 add_user_genus = False
            if add_user_genus:
                 priority_keywords.append(user_genus_equivalent)

    priority_keywords = list(filter(None, dict.fromkeys(priority_keywords))) 
    print(f"[WIKI_IMG] Priority keywords for image filename: {priority_keywords}")

    processed_search_terms = set() 

    while search_candidates:
        term_to_search = search_candidates.pop(0) 
        if not term_to_search or term_to_search in processed_search_terms:
            continue
        processed_search_terms.add(term_to_search)
        
        print(f"[WIKI_IMG] Trying Wikipedia search for term: '{term_to_search}'")
        try:
            wiki_page = None
            try:
                wiki_page = wikipedia.page(term_to_search, auto_suggest=True, redirect=True)
                print(f"[WIKI_IMG] Found page directly: '{wiki_page.title}' for '{term_to_search}'")
            except wikipedia.exceptions.PageError:
                print(f"[WIKI_IMG] Page not found directly for '{term_to_search}'. Trying wikipedia.search().")
                search_results = wikipedia.search(term_to_search, results=1)
                if not search_results:
                    print(f"[WIKI_IMG] No search results in Wikipedia for: '{term_to_search}'")
                    continue
                page_title_to_get = search_results[0]
                if page_title_to_get in processed_search_terms: 
                    print(f"[WIKI_IMG] Page '{page_title_to_get}' already processed, skipping.")
                    continue
                print(f"[WIKI_IMG] Found page via search: '{page_title_to_get}' for '{term_to_search}'")
                wiki_page = wikipedia.page(page_title_to_get, auto_suggest=False, redirect=True)
            
            if wiki_page and wiki_page.images:
                print(f"[WIKI_IMG] Page '{wiki_page.title}' images (up to 5): {wiki_page.images[:5]}")
                
                candidate_images_with_scores = []
                for img_url in wiki_page.images:
                    img_url_lower = img_url.lower()
                    if not any(ext in img_url_lower for ext in ['.png', '.jpg', '.jpeg']): 
                        continue
                    if any(keyword in img_url_lower for keyword in avoid_keywords_in_filename):
                        continue
                    
                    score = 0
                    for pk_word in priority_keywords:
                        if pk_word in img_url_lower:
                            score += 5 
                            filename_part = img_url_lower.split('/')[-1]
                            if filename_part.startswith(pk_word): 
                                score += 3
                            if pk_word == clean_scientific_name_for_filename and pk_word in filename_part : 
                                score +=5

                    if img_url_lower.endswith('.svg'): score -= 1 
                    candidate_images_with_scores.append({'url': img_url, 'score': score})
                
                if not candidate_images_with_scores:
                    print(f"[WIKI_IMG] No images passed filter for '{wiki_page.title}'")
                    continue

                sorted_images = sorted(candidate_images_with_scores, key=lambda x: x['score'], reverse=True)
                print(f"[WIKI_IMG] Sorted suitable images (top 3 with scores): {[{'url': i['url'][-50:], 'score': i['score']} for i in sorted_images[:3]]}")

                if sorted_images and sorted_images[0]['score'] > 0: 
                    best_image_url = sorted_images[0]['url']
                    if best_image_url.startswith("//"): best_image_url = "https:" + best_image_url
                    print(f"[WIKI_IMG] Best image found: {best_image_url} with score {sorted_images[0]['score']}")
                    return best_image_url
                else: 
                    print(f"[WIKI_IMG] No image found with positive score for '{wiki_page.title}'.")

            elif wiki_page:
                print(f"[WIKI_IMG] No images listed on Wikipedia page: '{wiki_page.title}'")
            
        except wikipedia.exceptions.DisambiguationError as e:
            print(f"[WIKI_IMG] Disambiguation for '{term_to_search}'. Options: {e.options[:3]}")
            if e.options:
                new_candidate = e.options[0]
                if new_candidate not in processed_search_terms and new_candidate not in search_candidates:
                    search_candidates.append(new_candidate)
                    print(f"[WIKI_IMG] Added disambiguation option '{new_candidate}' to search candidates.")
            continue 
        except wikipedia.exceptions.PageError:
             print(f"[WIKI_IMG] Wikipedia PageError (likely after search or disambiguation) for term: '{term_to_search}'")
             continue
        except Exception as e:
            print(f"[WIKI_IMG] Generic error for '{term_to_search}': {str(e)}")
            traceback.print_exc() 
            continue
            
    print(f"[WIKI_IMG] No suitable Wikipedia image URL after all attempts for initial candidates.")
    return None

# --- تابع جدید برای پیشنهاد نام علمی از NCBI ---
def get_best_ncbi_suggestion_flexible(common_name, max_ids_to_check=5):
    print(f"\n--- [NCBI Suggestion] Processing '{common_name}' ---")
    scientific_name_suggestion = None
    try:
        search_term = f"{common_name}[Common Name] OR {common_name}[Organism]"
        handle = Entrez.esearch(db="taxonomy", term=search_term, retmax=max_ids_to_check, sort="relevance")
        record = Entrez.read(handle)
        handle.close()
        id_list = record["IdList"]

        if not id_list:
            time.sleep(0.34) # NCBI rate limit compliance
            handle = Entrez.esearch(db="taxonomy", term=common_name, retmax=max_ids_to_check, sort="relevance")
            record = Entrez.read(handle)
            handle.close()
            id_list = record["IdList"]
            if not id_list:
                print(f"  [NCBI Suggestion] No TaxIDs found for '{common_name}' even after retry.")
                return None
        
        print(f"  [NCBI Suggestion] Found {len(id_list)} potential TaxIDs: {id_list}")
        ids_to_fetch = id_list[:max_ids_to_check] 
        if not ids_to_fetch: return None

        # NCBI rate limit: wait before next request (esummary)
        wait_time = 0.2 + (0.34 * len(ids_to_fetch)) 
        time.sleep(wait_time) 
        
        summary_handle = Entrez.esummary(db="taxonomy", id=",".join(ids_to_fetch), retmode="xml")
        summary_records = Entrez.read(summary_handle)
        summary_handle.close()
        
        candidates = { 
            "species": None, 
            "subspecies": None,
            "genus": None, 
            "family": None,
            "other_valid": None
        }
        
        for summary_record in summary_records:
            sci_name = summary_record.get("ScientificName")
            rank_ncbi = summary_record.get("Rank", "N/A").lower()
            # tax_id_processed = summary_record.get("Id") # برای لاگ کردن مفید است
            # print(f"    - Candidate TaxID {tax_id_processed}: '{sci_name}', Rank: {rank_ncbi}")

            if sci_name:
                if rank_ncbi == "species":
                    if not candidates["species"]:
                         candidates["species"] = sci_name
                elif rank_ncbi == "subspecies":
                    if not candidates["subspecies"]:
                        candidates["subspecies"] = sci_name
                elif rank_ncbi == "genus":
                     if not candidates["genus"]:
                         candidates["genus"] = sci_name
                elif rank_ncbi == "family":
                     if not candidates["family"]:
                         candidates["family"] = sci_name
                elif not candidates["other_valid"]: # اولین نام معتبر دیگر
                     candidates["other_valid"] = sci_name
        
        # تصمیم‌گیری نهایی بر اساس اولویت
        if candidates["species"]:
            scientific_name_suggestion = candidates["species"]
        elif candidates["subspecies"]:
            parts = candidates["subspecies"].split()
            if len(parts) >= 2: # استخراج نام گونه از زیرگونه
                scientific_name_suggestion = " ".join(parts[:2])
            else: # اگر فرمت نام زیرگونه استاندارد نبود، خود نام زیرگونه را برمی‌گردانیم
                scientific_name_suggestion = candidates["subspecies"] 
        elif candidates["genus"]:
            scientific_name_suggestion = candidates["genus"]
        elif candidates["family"]:
            scientific_name_suggestion = candidates["family"]
        elif candidates["other_valid"]:
             scientific_name_suggestion = candidates["other_valid"]
        # else:
            # print(f"  [NCBI Suggestion] No suitable scientific name found according to prioritization for '{common_name}'.")

    except Exception as e:
        print(f"  [NCBI Suggestion] Error querying NCBI Entrez for '{common_name}': {type(e).__name__} - {e}")
        traceback.print_exc() # برای دیباگ در لاگ‌های Vercel مفید است
        
    if scientific_name_suggestion:
        print(f"  [NCBI Suggestion] Final suggestion for '{common_name}': {scientific_name_suggestion}")
    else:
        print(f"  [NCBI Suggestion] No suggestion found for '{common_name}'.")
    return scientific_name_suggestion

# --- Endpoint جدید برای پیشنهاد نام علمی ---
@app.route('/api/suggest_scientific_name', methods=['GET', 'OPTIONS'])
def suggest_name_handler():
    common_headers = {'Access-Control-Allow-Origin': '*'}
    if request.method == 'OPTIONS':
        cors_headers = {
            **common_headers, 
            'Access-Control-Allow-Methods': 'GET, OPTIONS', 
            'Access-Control-Allow-Headers': 'Content-Type', 
            'Access-Control-Max-Age': '3600'
        }
        return ('', 204, cors_headers)

    common_name_query = request.args.get('common_name')
    if not common_name_query:
        return jsonify({"error": "پارامتر 'common_name' مورد نیاز است."}), 400, common_headers

    try:
        suggested_name = get_best_ncbi_suggestion_flexible(common_name_query)
        if suggested_name:
            return jsonify({
                "suggested_scientific_name": suggested_name, 
                "common_name_searched": common_name_query
            }), 200, common_headers
        else:
            return jsonify({
                "message": f"نام علمی برای '{common_name_query}' توسط NCBI پیشنهاد نشد.", 
                "common_name_searched": common_name_query
            }), 404, common_headers
    except Exception as e:
        print(f"Error in suggest_name_handler for '{common_name_query}': {str(e)}")
        traceback.print_exc()
        return jsonify({"error": "خطای داخلی سرور هنگام پردازش درخواست پیشنهاد نام."}), 500, common_headers

# --- Endpoint اصلی برای جستجوی طبقه‌بندی (با تغییرات جزئی در CORS برای سازگاری) ---
@app.route('/', defaults={'path': ''}, methods=['GET', 'POST', 'OPTIONS'])
@app.route('/<path:path>', methods=['GET', 'POST', 'OPTIONS'])
def main_handler(path=None):
    # اگر path برابر با 'api/suggest_scientific_name' بود، اینجا نباید اجرا شود
    # Flask ابتدا route های خاص تر را match می‌کند.
    # اما برای اطمینان، می‌توان یک شرط گذاشت یا مطمئن شد که Vercel routing درست کار می‌کند.
    # اگر endpoint های شما در یک فایل هستند، Flask خودش مدیریت می‌کند.
    # اگر `path` با `/api/suggest_scientific_name` شروع شود، این یعنی route بالا آن را نگرفته.
    # این اتفاق نباید بیفتد اگر endpoint ها درست تعریف شده باشند.

    common_headers = {'Access-Control-Allow-Origin': '*'}
    if request.method == 'OPTIONS':
        cors_headers = {
            **common_headers, 
            'Access-Control-Allow-Methods': 'GET, POST, OPTIONS', 
            'Access-Control-Allow-Headers': 'Content-Type, Authorization', 
            'Access-Control-Max-Age': '3600'
        }
        return ('', 204, cors_headers)

    species_name_query = ""
    if request.method == 'GET':
        species_name_query = request.args.get('name')
    elif request.method == 'POST':
        try:
            data_post = request.get_json()
            if data_post: species_name_query = data_post.get('name')
            else: return jsonify({"error": "درخواست نامعتبر: بدنه JSON خالی است یا قابل خواندن نیست."}), 400, common_headers
        except Exception as e_post:
            print(f"Error parsing JSON body: {str(e_post)}\n{traceback.format_exc()}")
            return jsonify({"error": f"خطا در پردازش درخواست: {str(e_post)}"}), 400, common_headers

    if not species_name_query:
        return jsonify({"error": "پارامتر 'name' (در آدرس یا بدنه JSON) مورد نیاز است."}), 400, common_headers

    params_gbif = {"name": species_name_query, "verbose": "true"}
    classification_data = {}
    gbif_error_message = None
    gbif_scientific_name = None

    try:
        api_response_gbif = requests.get(GBIF_API_URL_MATCH, params=params_gbif, timeout=10)
        api_response_gbif.raise_for_status()
        data_gbif = api_response_gbif.json()
        
        if not data_gbif or data_gbif.get("matchType") == "NONE" or data_gbif.get("confidence", 0) < 30:
            gbif_error_message = f"موجودی با نام '{species_name_query}' در GBIF پیدا نشد یا نتیجه با اطمینان کافی نبود."
            classification_data = {"searchedName": species_name_query, "matchType": data_gbif.get("matchType", "NONE"), "confidence": data_gbif.get("confidence")}
        else:
            gbif_scientific_name = data_gbif.get("scientificName")
            classification_data = {
                "searchedName": species_name_query,
                "scientificName": gbif_scientific_name,
                "kingdom": data_gbif.get("kingdom"), "phylum": data_gbif.get("phylum"), "class": data_gbif.get("class"),
                "order": data_gbif.get("order"), "family": data_gbif.get("family"), "genus": data_gbif.get("genus"),
                "species": data_gbif.get("species") if data_gbif.get("speciesKey") and data_gbif.get("species") else None,
                "usageKey": data_gbif.get("usageKey"), "confidence": data_gbif.get("confidence"),
                "matchType": data_gbif.get("matchType"), "status": data_gbif.get("status"), "rank": data_gbif.get("rank")
            }    
    except requests.exceptions.Timeout:
        gbif_error_message = "خطا: زمان پاسخگویی از سرور GBIF بیش از حد طول کشید."
        print(f"[GBIF_ERR] Timeout for: {species_name_query}")
    except requests.exceptions.HTTPError as http_err_gbif:
        gbif_error_message = f"خطا از سرور GBIF: {http_err_gbif}"
        try:
            gbif_error_details = api_response_gbif.json()
            gbif_error_message += f" - پیام GBIF: {gbif_error_details.get('message', api_response_gbif.text[:100])}"
        except: pass
        print(f"[GBIF_ERR] HTTPError: {gbif_error_message}")
    except requests.exceptions.RequestException as e_gbif_req:
        gbif_error_message = f"خطا در ارتباط با سرور GBIF: {str(e_gbif_req)}"
        print(f"[GBIF_ERR] RequestException: {str(e_gbif_req)}")
    except Exception as e_gbif_generic:
        gbif_error_message = "یک خطای پیش‌بینی نشده داخلی در ارتباط با GBIF رخ داده است."
        print(f"[GBIF_ERR] Generic Error: {str(e_gbif_generic)}")
        traceback.print_exc()

    wiki_image_url = get_wikipedia_image_url(species_name_query, gbif_scientific_name)
    
    if wiki_image_url:
        classification_data["imageUrl"] = wiki_image_url
    
    final_data = {k: v for k, v in classification_data.items() if v is not None}

    if gbif_error_message and not final_data.get("scientificName"): 
        if final_data.get("imageUrl"): # اگر GBIF خطا داد ولی تصویری پیدا شد (مثلا برای نام رایج)
             return jsonify({"message": gbif_error_message, **final_data}), 200, common_headers # با کد 200 چون تصویر داریم
        return jsonify({"message": gbif_error_message, "searchedName": species_name_query}), 404, common_headers
    
    return jsonify(final_data), 200, common_headers

# if __name__ == "__main__":
#    app.run(debug=True) # برای تست محلی، اگر لازم بود
