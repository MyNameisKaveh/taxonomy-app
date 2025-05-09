from flask import Flask, jsonify, request
import requests
import traceback
import wikipedia
# Entrez را دیرتر import می‌کنیم، بعد از پیکربندی
# from Bio import Entrez 
import time
import os 
import tempfile 

# تابع برای پیکربندی Entrez قبل از استفاده
def configure_entrez_for_readonly_env():
    try:
        from Bio import Entrez # Import Entrez اینجا
        from Bio.Entrez import Parser as EntrezParser

        # 1. تنظیم ایمیل (باید قبل از هر درخواست Entrez باشد)
        Entrez.email = "andolini1889@gmail.com" # ایمیل شما اینجا جایگذاری شد
        print(f"[CONFIG] Entrez.email set to: {Entrez.email}")

        # 2. تنظیم مسیر کش DTD
        custom_dtd_dir = os.path.join(tempfile.gettempdir(), "biopython_dtd_cache")
        if not os.path.exists(custom_dtd_dir):
            os.makedirs(custom_dtd_dir, mode=0o755, exist_ok=True) 
        
        EntrezParser.DataHandler.local_dtd_dir = custom_dtd_dir
        print(f"[CONFIG] Bio.Entrez DTD cache directory set to: {EntrezParser.DataHandler.local_dtd_dir}")
        
        return Entrez # برگرداندن ماژول Entrez پیکربندی شده

    except ImportError as ie:
        print(f"[CONFIG_FATAL] Could not import Bio.Entrez or Bio.Entrez.Parser: {ie}")
        raise 
    except Exception as e:
        print(f"[CONFIG_FATAL] Critical error setting Bio.Entrez DTD cache path: {type(e).__name__} - {e}")
        raise 

# Entrez را از طریق تابع پیکربندی دریافت می‌کنیم
Entrez = configure_entrez_for_readonly_env()


app = Flask(__name__)
GBIF_API_URL_MATCH = "https://api.gbif.org/v1/species/match"

# Entrez.email قبلا در configure_entrez_for_readonly_env() تنظیم شده است.

# --- تابع جستجوی تصویر از ویکی‌پدیا ---
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
        return None

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
    processed_search_terms = set() 
    while search_candidates:
        term_to_search = search_candidates.pop(0) 
        if not term_to_search or term_to_search in processed_search_terms:
            continue
        processed_search_terms.add(term_to_search)
        try:
            wiki_page = None
            try:
                wiki_page = wikipedia.page(term_to_search, auto_suggest=True, redirect=True)
            except wikipedia.exceptions.PageError:
                search_results = wikipedia.search(term_to_search, results=1)
                if not search_results: continue
                page_title_to_get = search_results[0]
                if page_title_to_get in processed_search_terms: continue
                wiki_page = wikipedia.page(page_title_to_get, auto_suggest=False, redirect=True)
            if wiki_page and wiki_page.images:
                candidate_images_with_scores = []
                for img_url in wiki_page.images:
                    img_url_lower = img_url.lower()
                    if not any(ext in img_url_lower for ext in ['.png', '.jpg', '.jpeg']): continue
                    if any(keyword in img_url_lower for keyword in avoid_keywords_in_filename): continue
                    score = 0
                    for pk_word in priority_keywords:
                        if pk_word in img_url_lower:
                            score += 5 
                            filename_part = img_url_lower.split('/')[-1]
                            if filename_part.startswith(pk_word): score += 3
                            if pk_word == clean_scientific_name_for_filename and pk_word in filename_part : score +=5
                    if img_url_lower.endswith('.svg'): score -= 1 
                    candidate_images_with_scores.append({'url': img_url, 'score': score})
                if not candidate_images_with_scores: continue
                sorted_images = sorted(candidate_images_with_scores, key=lambda x: x['score'], reverse=True)
                if sorted_images and sorted_images[0]['score'] > 0: 
                    best_image_url = sorted_images[0]['url']
                    if best_image_url.startswith("//"): best_image_url = "https:" + best_image_url
                    print(f"[WIKI_IMG] Best image for '{term_to_search}': {best_image_url} (Score: {sorted_images[0]['score']})")
                    return best_image_url
        except wikipedia.exceptions.DisambiguationError as e:
            if e.options:
                new_candidate = e.options[0]
                if new_candidate not in processed_search_terms and new_candidate not in search_candidates:
                    search_candidates.append(new_candidate)
            continue 
        except wikipedia.exceptions.PageError: continue
        except Exception as e:
            print(f"[WIKI_IMG_ERR] Generic for '{term_to_search}': {type(e).__name__} - {str(e)}")
            continue
    print(f"[WIKI_IMG] No image for user input: '{species_name_from_user}'.")
    return None

# --- تابع پیشنهاد نام علمی از NCBI ---
def get_best_ncbi_suggestion_flexible(common_name, max_ids_to_check=5):
    print(f"[NCBI_SUGGEST] Processing '{common_name}' (Email: {Entrez.email if Entrez else 'Email Not Set In Entrez Module'})")
    scientific_name_suggestion = None
    if not Entrez: 
        print("[NCBI_SUGGEST_ERR] Entrez module not available due to configuration error.")
        return None
    if not Entrez.email: # بررسی اضافه برای اطمینان از تنظیم بودن ایمیل
        print("[NCBI_SUGGEST_ERR] Entrez.email is not set. Please provide a valid email.")
        return None
        
    try:
        time.sleep(0.4) 
        search_term = f"{common_name}[Common Name] OR {common_name}[Organism]"
        handle = Entrez.esearch(db="taxonomy", term=search_term, retmax=max_ids_to_check, sort="relevance")
        record = Entrez.read(handle) 
        handle.close()
        id_list = record["IdList"]

        if not id_list:
            time.sleep(0.4)
            handle = Entrez.esearch(db="taxonomy", term=common_name, retmax=max_ids_to_check, sort="relevance")
            record = Entrez.read(handle)
            handle.close()
            id_list = record["IdList"]
            if not id_list:
                print(f"  [NCBI_SUGGEST] No TaxIDs for '{common_name}' after retry.")
                return None
        
        ids_to_fetch = id_list[:max_ids_to_check] 
        if not ids_to_fetch: return None

        time.sleep(0.4)
        summary_handle = Entrez.esummary(db="taxonomy", id=",".join(ids_to_fetch), retmode="xml")
        summary_records = Entrez.read(summary_handle)
        summary_handle.close()
        
        candidates = {"species": None, "subspecies": None, "genus": None, "family": None, "other_valid": None}
        for summary_record in summary_records:
            sci_name = summary_record.get("ScientificName")
            rank_ncbi = summary_record.get("Rank", "N/A").lower()
            if sci_name:
                if rank_ncbi == "species" and not candidates["species"]: candidates["species"] = sci_name
                elif rank_ncbi == "subspecies" and not candidates["subspecies"]: candidates["subspecies"] = sci_name
                elif rank_ncbi == "genus" and not candidates["genus"]: candidates["genus"] = sci_name
                elif rank_ncbi == "family" and not candidates["family"]: candidates["family"] = sci_name
                elif not candidates["other_valid"]: candidates["other_valid"] = sci_name
        
        if candidates["species"]: scientific_name_suggestion = candidates["species"]
        elif candidates["subspecies"]:
            parts = candidates["subspecies"].split()
            scientific_name_suggestion = " ".join(parts[:2]) if len(parts) >= 2 else candidates["subspecies"]
        elif candidates["genus"]: scientific_name_suggestion = candidates["genus"]
        elif candidates["family"]: scientific_name_suggestion = candidates["family"]
        elif candidates["other_valid"]: scientific_name_suggestion = candidates["other_valid"]

    except Exception as e:
        print(f"  [NCBI_SUGGEST_ERR] Querying NCBI for '{common_name}': {type(e).__name__} - {e}")
        if isinstance(e, OSError) and hasattr(e, 'errno') and e.errno == 30: # Read-only file system
            traceback.print_exc() 
        elif "HTTP Error 400" in str(e) or "Bad Request" in str(e): # خطای رایج اگر ایمیل مشکل داشته باشد
            print("  [NCBI_SUGGEST_WARN] Received HTTP 400 Bad Request. This might be due to an issue with the Entrez email or request rate.")
            traceback.print_exc() # برای دیدن جزئیات خطای HTTP
        
    if scientific_name_suggestion:
        print(f"  [NCBI_SUGGEST] Suggestion for '{common_name}': {scientific_name_suggestion}")
    else:
        print(f"  [NCBI_SUGGEST] No suggestion for '{common_name}'.")
    return scientific_name_suggestion

# --- Endpoint ها ---
@app.route('/api/suggest_scientific_name', methods=['GET', 'OPTIONS'])
def suggest_name_handler():
    common_headers = {'Access-Control-Allow-Origin': '*'}
    if request.method == 'OPTIONS':
        cors_headers = {**common_headers, 'Access-Control-Allow-Methods': 'GET, OPTIONS', 'Access-Control-Allow-Headers': 'Content-Type, X-Requested-With', 'Access-Control-Max-Age': '3600'}
        return ('', 204, cors_headers)
    common_name_query = request.args.get('common_name')
    if not common_name_query:
        return jsonify({"error": "پارامتر 'common_name' مورد نیاز است."}), 400, common_headers
    try:
        suggested_name = get_best_ncbi_suggestion_flexible(common_name_query)
        if suggested_name:
            return jsonify({"suggested_scientific_name": suggested_name, "common_name_searched": common_name_query}), 200, common_headers
        else:
            return jsonify({"message": f"نام علمی برای '{common_name_query}' توسط NCBI پیشنهاد نشد.", "common_name_searched": common_name_query}), 404, common_headers
    except Exception as e:
        print(f"[API_ERR] suggest_name_handler for '{common_name_query}': {type(e).__name__} - {str(e)}")
        traceback.print_exc()
        return jsonify({"error": "خطای داخلی سرور هنگام پردازش پیشنهاد نام."}), 500, common_headers

@app.route('/', defaults={'path': ''}, methods=['GET', 'POST', 'OPTIONS'])
@app.route('/<path:path>', methods=['GET', 'POST', 'OPTIONS'])
def main_handler(path=None):
    common_headers = {'Access-Control-Allow-Origin': '*'}
    if request.method == 'OPTIONS':
        cors_headers = {**common_headers, 'Access-Control-Allow-Methods': 'GET, POST, OPTIONS', 'Access-Control-Allow-Headers': 'Content-Type, Authorization, X-Requested-With', 'Access-Control-Max-Age': '3600'}
        return ('', 204, cors_headers)
    species_name_query = ""
    if request.method == 'GET': species_name_query = request.args.get('name')
    elif request.method == 'POST':
        try:
            data_post = request.get_json()
            if data_post: species_name_query = data_post.get('name')
            else: return jsonify({"error": "درخواست نامعتبر: بدنه JSON خالی."}), 400, common_headers
        except: return jsonify({"error": "خطا در خواندن بدنه JSON."}), 400, common_headers
    if not species_name_query: return jsonify({"error": "پارامتر 'name' نیاز است."}), 400, common_headers
    params_gbif = {"name": species_name_query, "verbose": "true"}
    classification_data = {}
    gbif_error_message = None
    gbif_scientific_name = None
    try:
        api_response_gbif = requests.get(GBIF_API_URL_MATCH, params=params_gbif, timeout=10)
        api_response_gbif.raise_for_status()
        data_gbif = api_response_gbif.json()
        if not data_gbif or data_gbif.get("matchType") == "NONE" or data_gbif.get("confidence", 0) < 30:
            gbif_error_message = f"موجودی '{species_name_query}' در GBIF پیدا نشد یا اطمینان کم بود."
            classification_data = {"searchedName": species_name_query, "matchType": data_gbif.get("matchType", "NONE"), "confidence": data_gbif.get("confidence")}
        else:
            gbif_scientific_name = data_gbif.get("scientificName")
            classification_data = {"searchedName": species_name_query, "scientificName": gbif_scientific_name, "kingdom": data_gbif.get("kingdom"), "phylum": data_gbif.get("phylum"), "class": data_gbif.get("class"), "order": data_gbif.get("order"), "family": data_gbif.get("family"), "genus": data_gbif.get("genus"), "species": data_gbif.get("species") if data_gbif.get("speciesKey") and data_gbif.get("species") else None, "usageKey": data_gbif.get("usageKey"), "confidence": data_gbif.get("confidence"), "matchType": data_gbif.get("matchType"), "status": data_gbif.get("status"), "rank": data_gbif.get("rank")}    
    except requests.exceptions.Timeout:
        gbif_error_message = "خطا: GBIF Timeout."
        print(f"[GBIF_ERR] Timeout: {species_name_query}")
    except requests.exceptions.HTTPError as http_err_gbif:
        gbif_error_message = f"خطا از GBIF: {http_err_gbif}"
        try:
            gbif_error_details = api_response_gbif.json()
            gbif_error_message += f" - پیام: {gbif_error_details.get('message', api_response_gbif.text[:100])}"
        except: pass
        print(f"[GBIF_ERR] HTTPError: {gbif_error_message}")
    except requests.exceptions.RequestException as e_gbif_req:
        gbif_error_message = f"خطا در ارتباط با GBIF: {str(e_gbif_req)}"
        print(f"[GBIF_ERR] RequestException: {str(e_gbif_req)}")
    except Exception as e_gbif_generic:
        gbif_error_message = "خطای داخلی GBIF."
        print(f"[GBIF_ERR] Generic Error: {str(e_gbif_generic)}")
        traceback.print_exc()
    wiki_image_url = get_wikipedia_image_url(species_name_query, gbif_scientific_name)
    if wiki_image_url: classification_data["imageUrl"] = wiki_image_url
    final_data = {k: v for k, v in classification_data.items() if v is not None}
    if gbif_error_message and not final_data.get("scientificName"): 
        if final_data.get("imageUrl"): return jsonify({"message": gbif_error_message, **final_data}), 200, common_headers
        return jsonify({"message": gbif_error_message, "searchedName": species_name_query}), 404, common_headers
    if not final_data and not gbif_error_message: return jsonify({"message": f"اطلاعاتی برای '{species_name_query}' یافت نشد."}), 404, common_headers
    return jsonify(final_data), 200, common_headers

# if __name__ == "__main__":
#    app.run(debug=True, port=5003)
