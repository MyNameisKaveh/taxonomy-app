from flask import Flask, jsonify, request
import requests
import traceback
import wikipedia # کتابخانه جدید برای کار با ویکی‌پدیا

app = Flask(__name__)
GBIF_API_URL_MATCH = "https://api.gbif.org/v1/species/match"

# تابع به‌روز شده برای گرفتن تصویر از ویکی‌پدیا
def get_wikipedia_image_url(species_name_from_user, scientific_name_from_gbif=None):
    search_candidates = []
    clean_scientific_name_for_filename = None

    if scientific_name_from_gbif:
        clean_scientific_name = scientific_name_from_gbif.split('(')[0].strip()
        if clean_scientific_name:
            search_candidates.append(clean_scientific_name)
            clean_scientific_name_for_filename = clean_scientific_name.lower().replace(" ", "_")
    
    user_name_for_filename = None
    if species_name_from_user:
        # اگر نام کاربر با نام علمی یکی بود، دوباره اضافه نکن
        if not (clean_scientific_name and species_name_from_user.lower() == clean_scientific_name.lower()):
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
        if "_" in clean_scientific_name_for_filename: # اضافه کردن بخش جنس به تنهایی
            priority_keywords.append(clean_scientific_name_for_filename.split("_")[0])
            
    if user_name_for_filename and user_name_for_filename not in priority_keywords:
        priority_keywords.append(user_name_for_filename)
        if "_" in user_name_for_filename: # اگر نام کاربر چند کلمه‌ای بود
             priority_keywords.append(user_name_for_filename.split("_")[0])


    priority_keywords = list(filter(None, dict.fromkeys(priority_keywords))) # حذف تکراری و None
    print(f"[WIKI_IMG] Priority keywords for image filename: {priority_keywords}")

    processed_search_terms = set() # برای جلوگیری از جستجوی تکراری یک عبارت

    while search_candidates:
        term_to_search = search_candidates.pop(0) # اولین کاندید رو بگیر و از لیست حذف کن
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
                if page_title_to_get in processed_search_terms: # اگر این صفحه قبلا پردازش شده
                    print(f"[WIKI_IMG] Page '{page_title_to_get}' already processed, skipping.")
                    continue
                print(f"[WIKI_IMG] Found page via search: '{page_title_to_get}' for '{term_to_search}'")
                wiki_page = wikipedia.page(page_title_to_get, auto_suggest=False, redirect=True)
            
            if wiki_page and wiki_page.images:
                print(f"[WIKI_IMG] Page '{wiki_page.title}' images (up to 5): {wiki_page.images[:5]}")
                
                candidate_images_with_scores = []
                for img_url in wiki_page.images:
                    img_url_lower = img_url.lower()
                    if not any(ext in img_url_lower for ext in ['.png', '.jpg', '.jpeg']): # .gif و .svg رو فعلا حذف کردم برای سادگی
                        continue
                    if any(keyword in img_url_lower for keyword in avoid_keywords_in_filename):
                        continue
                    
                    score = 0
                    # امتیاز بر اساس کلمات کلیدی اولویت‌دار
                    for pk_word in priority_keywords:
                        if pk_word in img_url_lower:
                            score += 5 # امتیاز بیشتر برای کلمه کلیدی دقیق
                            # اگر کلمه کلیدی در ابتدای نام فایل باشد (بعد از آخرین /)
                            filename_part = img_url_lower.split('/')[-1]
                            if filename_part.startswith(pk_word):
                                score += 3

                    # امتیاز کمتر برای svg اگر تصاویر دیگر هم هستند
                    if img_url_lower.endswith('.svg'):
                        score -= 1
                    
                    candidate_images_with_scores.append({'url': img_url, 'score': score})
                
                if not candidate_images_with_scores:
                    print(f"[WIKI_IMG] No images passed filter for '{wiki_page.title}'")
                    continue

                # مرتب‌سازی تصاویر بر اساس امتیاز (بیشترین امتیاز اول)
                sorted_images = sorted(candidate_images_with_scores, key=lambda x: x['score'], reverse=True)
                print(f"[WIKI_IMG] Sorted suitable images (top 3): {sorted_images[:3]}")

                best_image_url = sorted_images[0]['url']
                if best_image_url.startswith("//"): best_image_url = "https:" + best_image_url
                print(f"[WIKI_IMG] Best image found: {best_image_url} with score {sorted_images[0]['score']}")
                return best_image_url

            elif wiki_page:
                print(f"[WIKI_IMG] No images listed on Wikipedia page: '{wiki_page.title}'")
            
        except wikipedia.exceptions.DisambiguationError as e:
            print(f"[WIKI_IMG] Disambiguation for '{term_to_search}'. Options: {e.options[:3]}")
            if e.options:
                # اولین گزینه ابهام‌زدایی رو به انتهای لیست کاندیدها اضافه می‌کنیم اگر قبلا نبوده
                new_candidate = e.options[0]
                if new_candidate not in processed_search_terms and new_candidate not in search_candidates:
                    search_candidates.append(new_candidate)
                    print(f"[WIKI_IMG] Added disambiguation option '{new_candidate}' to search candidates.")
            continue # ادامه با کاندید بعدی یا گزینه ابهام‌زدایی
        except wikipedia.exceptions.PageError:
             print(f"[WIKI_IMG] Wikipedia PageError (likely after search or disambiguation) for term: '{term_to_search}'")
             continue
        except Exception as e:
            print(f"[WIKI_IMG] Generic error for '{term_to_search}': {str(e)}")
            traceback.print_exc() # چاپ کامل خطا برای دیباگ
            continue
            
    print(f"[WIKI_IMG] No suitable Wikipedia image URL after all attempts for initial candidates.")
    return None

# --- تابع main_handler و بقیه کدها مثل قبل باقی می‌ماند ---
# (همان کدی که در پیام قبلی برای main_handler و import ها دادم)

@app.route('/', defaults={'path': ''}, methods=['GET', 'POST', 'OPTIONS'])
@app.route('/<path:path>', methods=['GET', 'POST', 'OPTIONS'])
def main_handler(path=None):
    common_headers = {'Access-Control-Allow-Origin': '*'}

    if request.method == 'OPTIONS':
        cors_headers = {**common_headers, 'Access-Control-Allow-Methods': 'GET, POST, OPTIONS', 'Access-Control-Allow-Headers': 'Content-Type, Authorization', 'Access-Control-Max-Age': '3600'}
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
        if final_data.get("imageUrl"):
             return jsonify({"message": gbif_error_message, **final_data}), 200, common_headers
        return jsonify({"message": gbif_error_message, "searchedName": species_name_query}), 404, common_headers
    
    return jsonify(final_data), 200, common_headers
