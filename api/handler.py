from flask import Flask, jsonify, request
import requests

# تعریف اپلیکیشن Flask در سطح بالای ماژول
app = Flask(__name__)

# تعریف ثابت‌ها
GBIF_API_URL = "https://api.gbif.org/v1/species/match"

# تعریف تابع handler که به عنوان route اصلی عمل می‌کند
# این دکوراتورها app را به عنوان هدف می‌شناسند
@app.route('/', defaults={'path': ''}, methods=['GET', 'POST', 'OPTIONS'])
@app.route('/<path:path>', methods=['GET', 'POST', 'OPTIONS'])
def main_handler(path=None): # اسم تابع رو به main_handler تغییر دادم برای وضوح بیشتر، اختیاریه
    if request.method == 'OPTIONS':
        headers = {
            'Access-Control-Allow-Origin': '*',
            'Access-Control-Allow-Methods': 'GET, POST, OPTIONS',
            'Access-Control-Allow-Headers': 'Content-Type, Authorization', # Authorization رو هم اضافه کردم محض احتیاط
            'Access-Control-Max-Age': '3600'
        }
        return ('', 204, headers)

    species_name = ""
    if request.method == 'GET':
        species_name = request.args.get('name')
    elif request.method == 'POST':
        try:
            data = request.get_json()
            if data: # بررسی کنیم که data None نباشه
                species_name = data.get('name')
            else: # اگر بدنه JSON خالی بود یا قابل parse نبود
                return jsonify({"error": "Invalid or empty JSON body"}), 400, {'Access-Control-Allow-Origin': '*'}
        except Exception as e: # خطای کلی‌تر برای get_json
            return jsonify({"error": f"Error parsing JSON body: {str(e)}"}), 400, {'Access-Control-Allow-Origin': '*'}

    if not species_name:
        return jsonify({"error": "Parameter 'name' (in query string or JSON body) is required"}), 400, {'Access-Control-Allow-Origin': '*'}

    params = {
        "name": species_name,
        "verbose": "true"
    }

    try:
        api_response = requests.get(GBIF_API_URL, params=params)
        api_response.raise_for_status()
        data = api_response.json()

        if "usageKey" not in data or data.get("matchType") == "NONE" or not data.get("scientificName"):
            return jsonify({"error": f"Species '{species_name}' not found or not specific enough in GBIF."}), 404, {'Access-Control-Allow-Origin': '*'}

        classification = {
            "searchedName": species_name,
            "scientificName": data.get("scientificName"),
            "kingdom": data.get("kingdom"),
            "phylum": data.get("phylum"),
            "class": data.get("class"),
            "order": data.get("order"),
            "family": data.get("family"),
            "genus": data.get("genus"),
            "species": data.get("speciesKey") and data.get("species"),
            "confidence": data.get("confidence"),
            "matchType": data.get("matchType")
        }
        classification_cleaned = {k: v for k, v in classification.items() if v is not None}

        return jsonify(classification_cleaned), 200, {'Access-Control-Allow-Origin': '*'}

    except requests.exceptions.HTTPError as http_err:
        error_detail = f"HTTP error from GBIF: {http_err}"
        try:
            error_detail += f" - Response: {api_response.text}"
        except:
            pass
        return jsonify({"error": error_detail}), api_response.status_code if api_response else 500, {'Access-Control-Allow-Origin': '*'}
    except requests.exceptions.RequestException as e:
        return jsonify({"error": f"Error connecting to GBIF API: {str(e)}"}), 503, {'Access-Control-Allow-Origin': '*'}
    except Exception as e:
        # برای دیدن خطای دقیق‌تر در لاگ‌های Vercel
        import traceback
        print(f"Unexpected error: {str(e)}")
        print(traceback.format_exc())
        return jsonify({"error": f"An unexpected internal error occurred. Please check logs."}), 500, {'Access-Control-Allow-Origin': '*'}

# بخش if __name__ == "__main__": را برای Vercel کاملاً حذف کنید یا کامنت کنید.
# چون Vercel خودش app را پیدا و اجرا می‌کند.
