from flask import Flask, jsonify, request
import requests
import traceback
import wikipedia # کتابخانه جدید برای کار با ویکی‌پدیا

app = Flask(__name__)
GBIF_API_URL_MATCH = "https://api.gbif.org/v1/species/match"

def get_wikipedia_image_url(species_name_from_user, scientific_name_from_gbif=None):
    """
    سعی می‌کند URL تصویر اصلی گونه را از ویکی‌پدیای انگلیسی پیدا کند.
    """
    search_candidates = []
    if scientific_name_from_gbif:
        # نام علمی از GBIF (بدون نویسنده و سال) اولویت دارد
        clean_scientific_name = scientific_name_from_gbif.split('(')[0].strip()
        if clean_scientific_name:
            search_candidates.append(clean_scientific_name)
    
    # نام اولیه که کاربر وارد کرده
    if species_name_from_user:
        search_candidates.append(species_name_from_user)
    
    # اگر هیچ نامی برای جستجو نبود
    if not search_candidates:
        print("No search terms provided for Wikipedia image search.")
        return None

    print(f"[WIKI_IMG] Attempting to find image for candidates: {search_candidates}")
    wikipedia.set_lang("en") # زبان را یکبار در ابتدا تنظیم می‌کنیم

    for term_to_search in search_candidates:
        if not term_to_search:
            continue
        
        print(f"[WIKI_IMG] Trying Wikipedia search for term: '{term_to_search}'")
        try:
            # ابتدا سعی می‌کنیم مستقیماً صفحه را با نام بگیریم (auto_suggest=True به پیدا کردن کمک می‌کند)
            try:
                wiki_page = wikipedia.page(term_to_search, auto_suggest=True, redirect=True)
                print(f"[WIKI_IMG] Found page directly: '{wiki_page.title}' for term '{term_to_search}'")
            except wikipedia.exceptions.PageError:
                # اگر مستقیماً پیدا نشد، از search استفاده می‌کنیم
                print(f"[WIKI_IMG] Page not found directly for '{term_to_search}'. Using wikipedia.search().")
                search_results = wikipedia.search(term_to_search, results=1)
                if not search_results:
                    print(f"[WIKI_IMG] No search results in Wikipedia for: '{term_to_search}'")
                    continue # اگر جستجو هم نتیجه‌ای نداشت، سراغ عبارت بعدی برو
                
                page_title_to_get = search_results[0]
                print(f"[WIKI_IMG] Found page via search: '{page_title_to_get}' for term '{term_to_search}'")
                wiki_page = wikipedia.page(page_title_to_get, auto_suggest=False, redirect=True) # اینجا دیگه auto_suggest لازم نیست

            # حالا که صفحه را داریم، تصاویر را بررسی می‌کنیم
            if wiki_page.images:
                print(f"[WIKI_IMG] Images found on page '{wiki_page.title}': (showing up to 5) {wiki_page.images[:5]}")
                for img_url in wiki_page.images:
                    if img_url.lower().endswith(('.png', '.jpg', '.jpeg', '.gif', '.svg')):
                        # اطمینان از اینکه URL کامل است
                        if img_url.startswith("//"):
                            img_url = "https:" + img_url
                        # کتابخانه wikipedia معمولا URL های کامل برمی‌گردونه
                        print(f"[WIKI_IMG] Suitable image URL found: {img_url}")
                        return img_url
                print(f"[WIKI_IMG] No direct image URL with common extension found in page images for: '{wiki_page.title}'")
            else:
                print(f"[WIKI_IMG] No images listed on Wikipedia page: '{wiki_page.title}'")
            
        except wikipedia.exceptions.DisambiguationError as e:
            print(f"[WIKI_IMG] Wikipedia disambiguation error for '{term_to_search}'. Options (up to 3): {e.options[:3]}")
            if e.options:
                try:
                    # سعی می‌کنیم اولین گزینه از صفحه ابهام‌زدایی را استفاده کنیم
                    disambiguated_page_title = e.options[0]
                    print(f"[WIKI_IMG] Trying first disambiguation option: '{disambiguated_page_title}'")
                    wiki_page = wikipedia.page(disambiguated_page_title, auto_suggest=False, redirect=True)
                    if wiki_page.images:
                        print(f"[WIKI_IMG] Images found on disambiguated page '{wiki_page.title}': (up to 5) {wiki_page.images[:5]}")
                        for img_url in wiki_page.images:
                            if img_url.lower().endswith(('.png', '.jpg', '.jpeg', '.gif', '.svg')):
                                if img_url.startswith("//"): img_url = "https:" + img_url
                                print(f"[WIKI_IMG] Suitable image URL from disambiguated page: {img_url}")
                                return img_url
                        print(f"[WIKI_IMG] No direct image URL on disambiguated page images for: '{wiki_page.title}'")
                    else:
                        print(f"[WIKI_IMG] No images on disambiguated Wikipedia page: '{wiki_page.title}'")
                except Exception as de:
                    print(f"[WIKI_IMG] Error processing disambiguation option '{e.options[0]}': {de}")
            continue # اگر در ابهام‌زدایی موفق نبودیم، سراغ کاندید بعدی میریم
        except wikipedia.exceptions.PageError: # این خطا ممکنه بعد از search هم رخ بده اگر عنوان نامعتبر باشه
             print(f"[WIKI_IMG] Wikipedia PageError (likely after search) for term: '{term_to_search}'")
             continue
        except Exception as e:
            print(f"[WIKI_IMG] Generic error fetching Wikipedia image for '{term_to_search}': {str(e)}")
            # traceback.print_exc() # برای دیباگ بیشتر در لاگ Vercel می‌تونید اینو فعال کنید
            continue
            
    print(f"[WIKI_IMG] Could not find any suitable Wikipedia image URL after trying all candidates: {search_candidates}")
    return None


@app.route('/', defaults={'path': ''}, methods=['GET', 'POST', 'OPTIONS'])
@app.route('/<path:path>', methods=['GET', 'POST', 'OPTIONS'])
def main_handler(path=None):
    common_headers = {'Access-Control-Allow-Origin': '*'}

    if request.method == 'OPTIONS':
        # ... (کد CORS مثل قبل) ...
        cors_headers = {**common_headers, 'Access-Control-Allow-Methods': 'GET, POST, OPTIONS', 'Access-Control-Allow-Headers': 'Content-Type, Authorization', 'Access-Control-Max-Age': '3600'}
        return ('', 204, cors_headers)

    species_name_query = ""
    if request.method == 'GET':
        species_name_query = request.args.get('name')
    elif request.method == 'POST':
        try:
            data_post = request.get_json() # تغییر نام متغیر برای جلوگیری از تداخل با data از GBIF
            if data_post: species_name_query = data_post.get('name')
            else: return jsonify({"error": "درخواست نامعتبر: بدنه JSON خالی است یا قابل خواندن نیست."}), 400, common_headers
        except Exception as e_post: # تغییر نام متغیر خطا
            print(f"Error parsing JSON body: {str(e_post)}\n{traceback.format_exc()}")
            return jsonify({"error": f"خطا در پردازش درخواست: {str(e_post)}"}), 400, common_headers

    if not species_name_query:
        return jsonify({"error": "پارامتر 'name' (در آدرس یا بدنه JSON) مورد نیاز است."}), 400, common_headers

    params_gbif = {"name": species_name_query, "verbose": "true"} # تغییر نام متغیر
    classification_data = {}
    gbif_error_message = None
    gbif_scientific_name = None # برای ارسال به تابع تصویر

    try:
        api_response_gbif = requests.get(GBIF_API_URL_MATCH, params=params_gbif, timeout=10) # تغییر نام متغیر
        api_response_gbif.raise_for_status()
        data_gbif = api_response_gbif.json() # تغییر نام متغیر

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
    except requests.exceptions.HTTPError as http_err_gbif: # تغییر نام متغیر
        gbif_error_message = f"خطا از سرور GBIF: {http_err_gbif}"
        try:
            gbif_error_details = api_response_gbif.json() # تغییر نام متغیر
            gbif_error_message += f" - پیام GBIF: {gbif_error_details.get('message', api_response_gbif.text[:100])}"
        except: pass
        print(f"[GBIF_ERR] HTTPError: {gbif_error_message}")
    except requests.exceptions.RequestException as e_gbif_req: # تغییر نام متغیر
        gbif_error_message = f"خطا در ارتباط با سرور GBIF: {str(e_gbif_req)}"
        print(f"[GBIF_ERR] RequestException: {str(e_gbif_req)}")
    except Exception as e_gbif_generic: # تغییر نام متغیر
        gbif_error_message = "یک خطای پیش‌بینی نشده داخلی در ارتباط با GBIF رخ داده است."
        print(f"[GBIF_ERR] Generic Error: {str(e_gbif_generic)}")

    # تلاش برای گرفتن تصویر از ویکی‌پدیا
    wiki_image_url = get_wikipedia_image_url(species_name_query, gbif_scientific_name)
    
    if wiki_image_url:
        classification_data["imageUrl"] = wiki_image_url
    
    final_data = {k: v for k, v in classification_data.items() if v is not None}

    if gbif_error_message and not final_data.get("scientificName"): 
        # اگر خطای GBIF داشتیم و هیچ اطلاعاتی هم از GBIF نگرفتیم
        if final_data.get("imageUrl"): # اما تصویر از ویکی پیدا کردیم
             return jsonify({"message": gbif_error_message, **final_data}), 200, common_headers # با کد 200 چون تصویر داریم
        return jsonify({"message": gbif_error_message, "searchedName": species_name_query}), 404, common_headers
    
    return jsonify(final_data), 200, common_headers
