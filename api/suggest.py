 ```python
 # api/suggest.py
 from flask import Flask, jsonify, request
 from .utils import get_best_ncbi_suggestion_flexible # ایمپورت از فایل utils
 import os # برای خواندن متغیر محیطی اگر ایمیل را آنجا ذخیره کنیم

 app = Flask(__name__)

 # تنظیم ایمیل Entrez (یا از متغیر محیطی بخوانید)
 # Entrez.email = os.environ.get("ENTREZ_EMAIL", "your_default_email@example.com") 
 # فعلا همان ایمیل داخل utils استفاده می‌شود اگر آنجا تنظیم شده باشد.

 @app.route('/api/suggest_name', methods=['GET']) 
 def suggest_name_endpoint():
     """Endpoint فقط برای پیشنهاد نام علمی."""
     common_headers = {'Access-Control-Allow-Origin': '*'}
     query = request.args.get('query')
     lang = request.args.get('lang', 'en') 

     # افزودن هدرهای CORS به همه پاسخ‌ها (حتی خطاها)
     def make_response(data, status_code):
          return jsonify(data), status_code, common_headers

     if not query:
         return make_response({"error": "Query parameter 'query' is required"}, 400)
     if lang != 'en':
          return make_response({"error": "Currently only English queries are supported."}, 400)

     try:
         # فراخوانی تابع راهنما از utils.py
         suggestion = get_best_ncbi_suggestion_flexible(query)

         if suggestion:
             return make_response({"query": query, "scientific_name_suggestion": suggestion}, 200)
         else:
             return make_response({"query": query, "message": f"No scientific name suggestion found for '{query}'."}, 404)
     except Exception as e:
         # لاگ کردن خطای داخلی اگر در تابع راهنما رخ دهد
         print(f"Error in suggest_name_endpoint calling helper: {e}")
         import traceback
         traceback.print_exc()
         return make_response({"error": "An internal error occurred while getting suggestion."}, 500)

 # بخش if __name__ == "__main__": لازم نیست
 ```
