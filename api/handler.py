from flask import Flask, jsonify, request
import requests
import traceback
import wikipedia # کتابخانه جدید برای کار با ویکی‌پدیا

app = Flask(__name__)
GBIF_API_URL_MATCH = "https://api.gbif.org/v1/species/match"

# تنظیم زبان ویکی‌پدیا (اختیاری، اگر خواستید از ویکی‌پدیای فارسی استفاده کنید)
# wikipedia.set_lang("fa") # فعلا انگلیسی رو در نظر میگیریم چون تصاویر بیشتری داره

def get_wikipedia_image_url(species_name, scientific_name_from_gbif=None):
    """
    سعی می‌کند URL تصویر اصلی گونه را از ویکی‌پدیا پیدا کند.
    ابتدا با نام علمی دقیق از GBIF (اگر موجود باشد) و سپس با نام جستجو شده اولیه.
    """
    search_terms = []
    if scientific_name_from_gbif:
        search_terms.append(scientific_name_from_gbif.split('(')[0].strip()) # حذف بخش (Linnaeus, 1758)
    search_terms.append(species_name) # نام اولیه جستجو شده

    for term in search_terms:
        if not term:
            continue
        try:
            # تنظیم زبان برای هر جستجو (اگر خواستید بین زبانها سوییچ کنید)
            # wikipedia.set_lang("en") # یا "fa"
            
            # جستجوی صفحه با پیشنهاد (ممکنه عنوان دقیق نباشه)
            page_title_suggestion = wikipedia.suggest(term)
            if not page_title_suggestion: # اگر پیشنهادی نبود
                page_titles = wikipedia.search(term, results=1)
                if not page_titles:
                    continue # اگر جستجو هم نتیجه‌ای نداشت، سراغ عبارت بعدی برو
                page_title_suggestion = page_titles[0]
            
            wiki_page = wikipedia.page(page_title_suggestion, auto_suggest=False, redirect=True) # از پیشنهاد دقیق استفاده کن و ریدایرکت‌ها رو دنبال کن
            
            # گرفتن اولین تصویر از لیست تصاویر صفحه
            # تصاویر در ویکی‌پدیا معمولا با پسوندهای رایج هستن
            # ما به دنبال URL اصلی تصویر هستیم، نه thumbnail
            if wiki_page.images:
                for img_url in wiki_page.images:
                    if img_url.lower().endswith(('.png', '.jpg', '.jpeg', '.gif', '.svg')):
                        # گاهی اوقات URL ها در فرمت خاصی هستند، سعی می‌کنیم به URL قابل استفاده تبدیل کنیم
                        # این بخش ممکنه نیاز به بهبود داشته باشه
                        # مثال: commons.wikimedia.org/wiki/Special:FilePath/Lion_waiting_in_Nambia.jpg
                        # باید به .../Lion_waiting_in_Nambia.jpg تبدیل بشه
                        # فعلا ساده‌ترین حالت رو در نظر میگیریم
                        return img_url # اولین تصویر پیدا شده
            
        except wikipedia.exceptions.PageError:
            print(f"Wikipedia page not found for: {term}")
            continue # اگر صفحه‌ای پیدا نشد، سراغ عبارت بعدی برو
        except wikipedia.exceptions.DisambiguationError as e:
            print(f"Wikipedia disambiguation error for {term}: {e.options[:3]}") # صفحه ابهام‌زدایی
            # می‌تونیم سعی کنیم اولین گزینه از صفحه ابهام‌زدایی رو انتخاب کنیم
            if e.options:
                try:
                    disambiguated_page = wikipedia.page(e.options[0], auto_suggest=False, redirect=True)
                    if disambiguated_page.images:
                         for img_url in disambiguated_page.images:
                            if img_url.lower().endswith(('.png', '.jpg', '.jpeg', '.gif', '.svg')):
                                return img_url
                except Exception as de:
                    print(f"Error getting image from disambiguated page: {de}")
            continue
        except Exception as e:
            print(f"Error fetching Wikipedia image for {term}: {str(e)}\n{traceback.format_exc()}")
            continue # در صورت بروز خطای دیگر، سراغ عبارت بعدی برو
    return None


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
        # ... (کد قبلی برای گرفتن نام از POST) ...
        try:
            data = request.get_json()
            if data: species_name_query = data.get('name')
            else: return jsonify({"error": "درخواست نامعتبر: بدنه JSON خالی است یا قابل خواندن نیست."}), 400, common_headers
        except Exception as e:
            print(f"Error parsing JSON body: {str(e)}\n{traceback.format_exc()}")
            return jsonify({"error": f"خطا در پردازش درخواست: {str(e)}"}), 400, common_headers


    if not species_name_query:
        return jsonify({"error": "پارامتر 'name' (در آدرس یا بدنه JSON) مورد نیاز است."}), 400, common_headers

    params = {"name": species_name_query, "verbose": "true"}

    classification_data = {}
    gbif_error_message = None

    try:
        api_response = requests.get(GBIF_API_URL_MATCH, params=params, timeout=10)
        api_response.raise_for_status()
        data = api_response.json()

        if not data or data.get("matchType") == "NONE" or data.get("confidence", 0) < 30:
            gbif_error_message = f"موجودی با نام '{species_name_query}' در GBIF پیدا نشد یا نتیجه با اطمینان کافی نبود."
            classification_data = {"searchedName": species_name_query, "matchType": data.get("matchType", "NONE"), "confidence": data.get("confidence")}
        else:
            classification_data = {
                "searchedName": species_name_query,
                "scientificName": data.get("scientificName"),
                "kingdom": data.get("kingdom"), "phylum": data.get("phylum"), "class": data.get("class"),
                "order": data.get("order"), "family": data.get("family"), "genus": data.get("genus"),
                "species": data.get("species") if data.get("speciesKey") and data.get("species") else None,
                "usageKey": data.get("usageKey"), "confidence": data.get("confidence"),
                "matchType": data.get("matchType"), "status": data.get("status"), "rank": data.get("rank")
            }
    
    except requests.exceptions.Timeout:
        gbif_error_message = "خطا: زمان پاسخگویی از سرور GBIF بیش از حد طول کشید."
        print(f"Timeout connecting to GBIF for: {species_name_query}")
    except requests.exceptions.HTTPError as http_err:
        # ... (کد قبلی برای مدیریت خطای HTTP از GBIF) ...
        gbif_error_message = f"خطا از سرور GBIF: {http_err}"
        # ... (بقیه کد)
    except requests.exceptions.RequestException as e:
        gbif_error_message = f"خطا در ارتباط با سرور GBIF: {str(e)}"
        print(f"Network error connecting to GBIF: {str(e)}\n{traceback.format_exc()}")
    except Exception as e:
        gbif_error_message = "یک خطای پیش‌بینی نشده داخلی در ارتباط با GBIF رخ داده است."
        print(f"Unexpected GBIF error: {str(e)}\n{traceback.format_exc()}")

    # تلاش برای گرفتن تصویر از ویکی‌پدیا حتی اگر GBIF خطا داده باشه (با نام اولیه)
    # یا اگر GBIF نتیجه داده، با نام علمی دقیق‌تر
    wiki_image_url = None
    if classification_data.get("scientificName"):
        wiki_image_url = get_wikipedia_image_url(species_name_query, classification_data["scientificName"])
    else: # اگر نام علمی از GBIF نداشتیم، با نام اولیه جستجو کن
        wiki_image_url = get_wikipedia_image_url(species_name_query)
    
    if wiki_image_url:
        classification_data["imageUrl"] = wiki_image_url
    
    # برگرداندن نتیجه نهایی
    if gbif_error_message and not classification_data.get("scientificName"): # اگر خطای GBIF داشتیم و هیچ اطلاعاتی هم نگرفتیم
        # اگر imageUrl داشتیم، اون رو با پیام خطا برمیگردونیم
        if wiki_image_url:
             return jsonify({"message": gbif_error_message, **classification_data}), 404, common_headers
        return jsonify({"message": gbif_error_message, "searchedName": species_name_query}), 404, common_headers
    
    # پاک کردن کلیدهایی که مقدار ندارند از classification_data نهایی
    final_data = {k: v for k, v in classification_data.items() if v is not None}
    return jsonify(final_data), 200, common_headers
