from flask import Flask, jsonify, request
import requests
import traceback # برای لاگ کردن خطاهای دقیق‌تر

app = Flask(__name__)
GBIF_API_URL_MATCH = "https://api.gbif.org/v1/species/match"
# GBIF_API_URL_SPECIES = "https://api.gbif.org/v1/species/" # برای گرفتن اطلاعات بیشتر و تصویر در مرحله بعد

@app.route('/', defaults={'path': ''}, methods=['GET', 'POST', 'OPTIONS'])
@app.route('/<path:path>', methods=['GET', 'POST', 'OPTIONS'])
def main_handler(path=None):
    common_headers = {'Access-Control-Allow-Origin': '*'} # هدر مشترک

    if request.method == 'OPTIONS':
        cors_headers = {
            **common_headers, # اضافه کردن هدر مشترک
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
            data = request.get_json()
            if data:
                species_name_query = data.get('name')
            else:
                return jsonify({"error": "درخواست نامعتبر: بدنه JSON خالی است یا قابل خواندن نیست."}), 400, common_headers
        except Exception as e:
            print(f"Error parsing JSON body: {str(e)}\n{traceback.format_exc()}")
            return jsonify({"error": f"خطا در پردازش درخواست: {str(e)}"}), 400, common_headers

    if not species_name_query:
        return jsonify({"error": "پارامتر 'name' (در آدرس یا بدنه JSON) مورد نیاز است."}), 400, common_headers

    # پارامترهای پایه برای GBIF
    # verbose=true برای اطلاعات کامل‌تر از جمله طبقه‌بندی‌های بالاتر
    # kingdom, phylum, class, order, family, genus, species رو معمولا برمیگردونه اگر پیدا کنه
    params = {
        "name": species_name_query,
        "verbose": "true"
    }

    try:
        api_response = requests.get(GBIF_API_URL_MATCH, params=params, timeout=10) # اضافه کردن timeout
        api_response.raise_for_status() # اگر خطای HTTP داشت (4xx, 5xx), exception ایجاد میکنه
        data = api_response.json()

        # بررسی اولیه برای اینکه آیا اصلا نتیجه‌ای هست یا نه
        if not data or data.get("matchType") == "NONE" or data.get("confidence", 0) < 30: # حداقل confidence رو 30 در نظر میگیریم
             # اگر confidence خیلی پایین بود یا matchType نبود، یعنی نتیجه مناسبی پیدا نشده
            return jsonify({
                "message": f"موجودی با نام '{species_name_query}' در GBIF پیدا نشد یا نتیجه با اطمینان کافی نبود.",
                "searchedName": species_name_query,
                "matchType": data.get("matchType", "NONE"),
                "confidence": data.get("confidence")
            }), 404, common_headers


        # استخراج اطلاعات طبقه‌بندی از نتیجه اصلی (بهترین تطابق)
        # GBIF ممکنه فیلدهای speciesKey و species رو فقط برای تطابق دقیق با گونه برگردونه
        # برای رتبه‌های بالاتر مثل جنس، این فیلدها ممکنه نباشن
        classification = {
            "searchedName": species_name_query,
            "scientificName": data.get("scientificName"),
            "kingdom": data.get("kingdom"),
            "phylum": data.get("phylum"),
            "class": data.get("class"),
            "order": data.get("order"),
            "family": data.get("family"),
            "genus": data.get("genus"),
            # اگر speciesKey وجود داره و species هم مقدار داره، اون رو برمیگردونیم
            "species": data.get("species") if data.get("speciesKey") and data.get("species") else None,
            "usageKey": data.get("usageKey"), # usageKey برای گرفتن اطلاعات بیشتر (مثلا تصویر) لازمه
            "confidence": data.get("confidence"),
            "matchType": data.get("matchType"),
            "status": data.get("status"), # وضعیت تاکسونومیکی (ACCEPTED, DOUBTFUL, SYNONYM)
            "rank": data.get("rank") # رتبه تاکسونومیکی نتیجه (SPECIES, GENUS, FAMILY, ...)
        }
        
        classification_cleaned = {k: v for k, v in classification.items() if v is not None}

        return jsonify(classification_cleaned), 200, common_headers

    except requests.exceptions.Timeout:
        print(f"Timeout connecting to GBIF for: {species_name_query}")
        return jsonify({"error": "خطا: زمان پاسخگویی از سرور GBIF بیش از حد طول کشید."}), 504, common_headers # Gateway Timeout
    except requests.exceptions.HTTPError as http_err:
        error_message = f"خطا از سرور GBIF: {http_err}"
        try:
            gbif_error_data = api_response.json()
            if isinstance(gbif_error_data, dict) and "message" in gbif_error_data:
                error_message += f" - پیام GBIF: {gbif_error_data['message']}"
            else:
                 error_message += f" - پاسخ GBIF: {api_response.text[:200]}" # فقط بخشی از پاسخ برای جلوگیری از لاگ طولانی
        except:
            error_message += f" - (عدم توانایی در خواندن پیام خطای GBIF: {api_response.status_code})"
        print(error_message)
        return jsonify({"error": error_message}), api_response.status_code if api_response else 500, common_headers
    except requests.exceptions.RequestException as e:
        print(f"Network error connecting to GBIF: {str(e)}\n{traceback.format_exc()}")
        return jsonify({"error": f"خطا در ارتباط با سرور GBIF: {str(e)}"}), 503, common_headers
    except Exception as e:
        print(f"Unexpected error in handler: {str(e)}\n{traceback.format_exc()}")
        return jsonify({"error": "یک خطای پیش‌بینی نشده داخلی رخ داده است. لطفاً با پشتیبانی تماس بگیرید یا لاگ‌ها را بررسی کنید."}), 500, common_headers
