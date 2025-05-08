# api/handler.py

from flask import Flask, jsonify, request
import requests
import traceback
import wikipedia 
from .utils import get_best_ncbi_suggestion_flexible # ایمپورت از فایل utils

app = Flask(__name__)
GBIF_API_URL_MATCH = "https://api.gbif.org/v1/species/match"

# =============================================
# تابع کمکی برای گرفتن تصویر از ویکی‌پدیا
# =============================================
def get_wikipedia_image_url(species_name_from_user, scientific_name_from_gbif=None):
    # ... (کد کامل و صحیح تابع get_wikipedia_image_url که قبلا داشتیم و تصویر درست رو برمی‌گردوند) ...
    # ... (مطمئن شو کد کامل این تابع اینجا هست) ...
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
        if not term_to_search or term_to_search in processed_search_terms: continue
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
                if not search_results: print(f"[WIKI_IMG] No search results in Wikipedia for: '{term_to_search}'"); continue
                page_title_to_get = search_results[0]
                if page_title_to_get in processed_search_terms: print(f"[WIKI_IMG] Page '{page_title_to_get}' already processed, skipping."); continue
                print(f"[WIKI_IMG] Found page via search: '{page_title_to_get}' for '{term_to_search}'")
                wiki_page = wikipedia.page(page_title_to_get, auto_suggest=False, redirect=True)            
            if wiki_page and wiki_page.images:
                print(f"[WIKI_IMG] Page '{wiki_page.title}' images (up to 5): {wiki_page.images[:5]}")                
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
                if not candidate_images_with_scores: print(f"[WIKI_IMG] No images passed filter for '{wiki_page.title}'"); continue
                sorted_images = sorted(candidate_images_with_scores, key=lambda x: x['score'], reverse=True)
                print(f"[WIKI_IMG] Sorted suitable images (top 3 with scores): {[{'url': i['url'][-50:], 'score': i['score']} for i in sorted_images[:3]]}")
                if sorted_images and sorted_images[0]['score'] > 0:
                    best_image_url = sorted_images[0]['url']
                    if best_image_url.startswith("//"): best_image_url = "https:" + best_image_url
                    print(f"[WIKI_IMG] Best image found: {best_image_url} with score {sorted_images[0]['score']}")
                    return best_image_url
                else: print(f"[WIKI_IMG] No image found with positive score for '{wiki_page.title}'.")
            elif wiki_page: print(f"[WIKI_IMG] No images listed on Wikipedia page: '{wiki_page.title}'")            
        except wikipedia.exceptions.DisambiguationError as e:
            print(f"[WIKI_IMG] Disambiguation for '{term_to_search}'. Options: {e.options[:3]}")
            if e.options:
                new_candidate = e.options[0]
                if new_candidate not in processed_search_terms and new_candidate not in search_candidates:
                    search_candidates.append(new_candidate)
                    print(f"[WIKI_IMG] Added disambiguation option '{new_candidate}' to search candidates.")
            continue 
        except wikipedia.exceptions.PageError: print(f"[WIKI_IMG] Wikipedia PageError (likely after search or disambiguation) for term: '{term_to_search}'"); continue
        except Exception as e: print(f"[WIKI_IMG] Generic error for '{term_to_search}': {str(e)}"); traceback.print_exc(); continue            
    print(f"[WIKI_IMG] No suitable Wikipedia image URL after all attempts for initial candidates.")
    return None

# =============================================
# Endpoint جدید برای راهنمای نام علمی
# =============================================
@app.route('/api/suggest_name', methods=['GET'])
def suggest_name_endpoint():
    """Endpoint برای پیشنهاد نام علمی بر اساس نام رایج."""
    common_headers = {'Access-Control-Allow-Origin': '*'}
    query = request.args.get('query')
    lang = request.args.get('lang', 'en') # فعلا فقط انگلیسی 'en' پشتیبانی می‌شود

    if not query:
        return jsonify({"error": "Query parameter 'query' is required"}), 400, common_headers
    if lang != 'en':
         return jsonify({"error": "Currently only English queries are supported."}), 400, common_headers

    # فراخوانی تابع راهنما از utils.py
    suggestion = get_best_ncbi_suggestion_flexible(query)

    if suggestion:
        return jsonify({"query": query, "scientific_name_suggestion": suggestion}), 200, common_headers
    else:
        # اگر در utils.py برای خطاها None برمیگردد، اینجا 404 میدهیم
        return jsonify({"query": query, "message": f"No scientific name suggestion found for '{query}'."}), 404, common_headers


# =============================================
# Endpoint اصلی برای جستجوی طبقه‌بندی
# =============================================
@app.route('/api/handler', methods=['GET', 'POST', 'OPTIONS']) # <-- تغییر مسیر اصلی به /api/handler
@app.route('/', defaults={'path': ''}, methods=['GET', 'POST', 'OPTIONS']) # <-- route ریشه هم به همین تابع بره
@app.route('/<path:path>', methods=['GET', 'POST', 'OPTIONS']) # <-- بقیه مسیرها هم (برای سازگاری)
def main_handler(path=None):
    """Endpoint اصلی برای جستجوی طبقه‌بندی و تصویر."""
    common_headers = {'Access-Control-Allow-Origin': '*'}
    if request.method == 'OPTIONS':
        cors_headers = {**common_headers, 'Access-Control-Allow-Methods': 'GET, POST, OPTIONS', 'Access-Control-Allow-Headers': 'Content-Type, Authorization', 'Access-Control-Max-Age': '3600'}
        return ('', 204, cors_headers)

    # --- گرفتن نام جستجو از query param یا body ---
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
         # اگر نامی نرسیده، شاید درخواست به ریشه بوده، یک صفحه پیش‌فرض یا راهنما نشان دهیم؟
         # فعلا همان خطای قبلی را میدهیم.
        return jsonify({"error": "پارامتر 'name' (در آدرس یا بدنه JSON) مورد نیاز است."}), 400, common_headers

    # --- جستجو در GBIF ---
    params_gbif = {"name": species_name_query, "verbose": "true"}
    classification_data = {"searchedName": species_name_query}
    gbif_error_message = None
    gbif_scientific_name = None
    try:
        api_response_gbif = requests.get(GBIF_API_URL_MATCH, params=params_gbif, timeout=10)
        api_response_gbif.raise_for_status()
        data_gbif = api_response_gbif.json()
        if not data_gbif or data_gbif.get("matchType") == "NONE" or data_gbif.get("confidence", 0) < 30:
            gbif_error_message = f"موجودی با نام '{species_name_query}' در GBIF پیدا نشد یا نتیجه با اطمینان کافی نبود."
            classification_data["matchType"] = data_gbif.get("matchType", "NONE")
            classification_data["confidence"] = data_gbif.get("confidence")
        else:
            gbif_scientific_name = data_gbif.get("scientificName")
            classification_data.update({
                "scientificName": gbif_scientific_name, "kingdom": data_gbif.get("kingdom"), 
                "phylum": data_gbif.get("phylum"), "class": data_gbif.get("class"),
                "order": data_gbif.get("order"), "family": data_gbif.get("family"), 
                "genus": data_gbif.get("genus"), 
                "species": data_gbif.get("species") if data_gbif.get("speciesKey") and data_gbif.get("species") else None,
                "usageKey": data_gbif.get("usageKey"), "confidence": data_gbif.get("confidence"),
                "matchType": data_gbif.get("matchType"), "status": data_gbif.get("status"), 
                "rank": data_gbif.get("rank")
            })   
    except requests.exceptions.Timeout:
        gbif_error_message = "خطا: زمان پاسخگویی از سرور GBIF بیش از حد طول کشید."; print(f"[GBIF_ERR] Timeout for: {species_name_query}")
    except requests.exceptions.HTTPError as http_err_gbif:
        gbif_error_message = f"خطا از سرور GBIF: {http_err_gbif}"; print(f"[GBIF_ERR] HTTPError: {gbif_error_message}")
    except requests.exceptions.RequestException as e_gbif_req:
        gbif_error_message = f"خطا در ارتباط با سرور GBIF: {str(e_gbif_req)}"; print(f"[GBIF_ERR] RequestException: {str(e_gbif_req)}")
    except Exception as e_gbif_generic:
        gbif_error_message = "یک خطای پیش‌بینی نشده داخلی در ارتباط با GBIF رخ داده است."; print(f"[GBIF_ERR] Generic Error: {str(e_gbif_generic)}"); traceback.print_exc()

    # --- جستجو برای تصویر در ویکی‌پدیا ---
    wiki_image_url = get_wikipedia_image_url(species_name_query, gbif_scientific_name)
    if wiki_image_url:
        classification_data["imageUrl"] = wiki_image_url
    
    # --- آماده‌سازی و بازگرداندن پاسخ نهایی ---
    final_data = {k: v for k, v in classification_data.items() if v is not None}
    if gbif_error_message and not final_data.get("scientificName"): 
        if final_data.get("imageUrl"):
             final_data["gbifLookupMessage"] = gbif_error_message 
             return jsonify(final_data), 200, common_headers
        return jsonify({"message": gbif_error_message, "searchedName": species_name_query}), 404, common_headers
    
    return jsonify(final_data), 200, common_headers
