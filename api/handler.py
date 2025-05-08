from flask import Flask, jsonify, request
import requests

app = Flask(__name__)

GBIF_API_URL = "https://api.gbif.org/v1/species/match"

@app.route('/', defaults={'path': ''}, methods=['GET', 'POST']) # Allow both GET and POST
@app.route('/<path:path>', methods=['GET', 'POST']) # Allow both GET and POST for any subpath
def handler(path=None): # Add path=None as default
    # برای اینکه Vercel تشخیص بده این یک تابع است و بتواند با هر روشی آن را فراخوانی کند
    # و هر مسیری داخل /api/handler را بگیرد.

    if request.method == 'OPTIONS': # Handle CORS preflight requests
        headers = {
            'Access-Control-Allow-Origin': '*',
            'Access-Control-Allow-Methods': 'GET, POST, OPTIONS',
            'Access-Control-Allow-Headers': 'Content-Type',
            'Access-Control-Max-Age': '3600'
        }
        return ('', 204, headers)

    # دریافت نام گونه از query parameter (برای تست با مرورگر) یا از JSON body (برای فرانت‌اند)
    species_name = ""
    if request.method == 'GET':
        species_name = request.args.get('name')
    elif request.method == 'POST':
        try:
            data = request.get_json()
            species_name = data.get('name')
        except:
            return jsonify({"error": "Invalid JSON body"}), 400


    if not species_name:
        return jsonify({"error": "Parameter 'name' (in query string or JSON body) is required"}), 400

    params = {
        "name": species_name,
        "verbose": "true" # برای دریافت اطلاعات بیشتر از جمله طبقه‌بندی بالاتر
        # "kingdom": "Animalia" # میتونید برای شروع به یک فرمانرو خاص محدود کنید
    }

    try:
        api_response = requests.get(GBIF_API_URL, params=params)
        api_response.raise_for_status() # اگر خطای HTTP داشت، exception ایجاد میکنه
        data = api_response.json()

        # بررسی اینکه آیا گونه‌ای پیدا شده یا نه
        # GBIF ممکنه یک نتیجه "matchType: NONE" برگردونه یا فیلدهای اصلی خالی باشن
        if "usageKey" not in data or data.get("matchType") == "NONE" or not data.get("scientificName"):
            return jsonify({"error": f"Species '{species_name}' not found or not specific enough in GBIF."}), 404

        classification = {
            "searchedName": species_name,
            "scientificName": data.get("scientificName"),
            "kingdom": data.get("kingdom"),
            "phylum": data.get("phylum"),
            "class": data.get("class"),
            "order": data.get("order"),
            "family": data.get("family"),
            "genus": data.get("genus"),
            "species": data.get("speciesKey") and data.get("species"), # گونه ممکنه در برخی سطوح بالاتر نباشه
            "confidence": data.get("confidence"),
            "matchType": data.get("matchType")
        }
        # حذف کلیدهایی که مقدار ندارند
        classification_cleaned = {k: v for k, v in classification.items() if v is not None}

        # اضافه کردن هدرهای CORS به پاسخ اصلی
        response_headers = {
            'Access-Control-Allow-Origin': '*'
        }
        return (jsonify(classification_cleaned), 200, response_headers)

    except requests.exceptions.HTTPError as http_err:
        # برای خطاهای HTTP از خود GBIF (مثل 404 اگر API به شکل دیگری خطا بدهد)
        error_detail = f"HTTP error from GBIF: {http_err}"
        try: # سعی در خواندن متن خطای GBIF
            error_detail += f" - Response: {api_response.text}"
        except:
            pass
        return jsonify({"error": error_detail}), api_response.status_code
    except requests.exceptions.RequestException as e:
        # برای خطاهای شبکه در ارتباط با GBIF
        return jsonify({"error": f"Error connecting to GBIF API: {str(e)}"}), 503 # Service Unavailable
    except Exception as e:
        # برای خطاهای پیش‌بینی نشده دیگر در کد ما
        return jsonify({"error": f"An unexpected error occurred: {str(e)}"}), 500

# این بخش برای اجرای محلی است و روی Vercel استفاده نمی‌شود،
# اما برای اینکه ساختار یک اپ Flask کامل باشه، نگهش می‌داریم.
# if __name__ == "__main__":
#     app.run(port=5000) # پورت را برای Vercel مشخص نمی‌کنیم
