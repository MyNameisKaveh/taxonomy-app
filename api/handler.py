# --- Imports ---
from flask import Flask, jsonify, request
import requests
import traceback
import wikipedia
from Bio import Entrez # برای NCBI
import time

# --- Flask App Initialization ---
app = Flask(__name__)

# --- Configuration ---
GBIF_API_URL_MATCH = "https://api.gbif.org/v1/species/match"
# !!! تنظیم ایمیل برای NCBI Entrez !!!
Entrez.email = "andolini1889@gmail.com" 

# --- Helper Function: NCBI Suggestion ---
def get_best_ncbi_suggestion_flexible(common_name_en, max_ids_to_check=5):
    """
    با استفاده از NCBI Entrez، بهترین نام علمی ممکن را برای یک نام رایج انگلیسی پیدا می‌کند.
    """
    print(f"\n--- Processing '{common_name_en}' with NCBI Entrez (Flexible Ranks) ---")
    scientific_name_suggestion = None
    try:
        search_term = common_name_en 
        print(f"  Searching NCBI Taxonomy with term: '{search_term}'")
        
        handle = Entrez.esearch(db="taxonomy", term=search_term, retmax=max_ids_to_check, sort="relevance")
        record = Entrez.read(handle)
        handle.close()
        id_list = record["IdList"]

        if not id_list:
            print(f"  No TaxIDs found for '{common_name_en}'.")
            return None
        
        print(f"  Found {len(id_list)} potential TaxIDs: {id_list}")
        ids_to_fetch = id_list[:max_ids_to_check] 
        if not ids_to_fetch: return None

        wait_time = 0.2 + (0.34 * len(ids_to_fetch)) 
        time.sleep(wait_time) 
        
        summary_handle = Entrez.esummary(db="taxonomy", id=",".join(ids_to_fetch), retmode="xml")
        summary_records = Entrez.read(summary_handle)
        summary_handle.close()
        
        candidates = {"species": None, "subspecies": None, "genus": None, "family": None, "other_valid": None}
        
        for summary_record in summary_records:
            sci_name = summary_record.get("ScientificName")
            rank_ncbi = summary_record.get("Rank", "N/A").lower()
            tax_id_processed = summary_record.get("Id") 
            # print(f"    - Candidate TaxID {tax_id_processed}: '{sci_name}', Rank: {rank_ncbi}") 
            if sci_name:
                if rank_ncbi == "species" and not candidates["species"]: candidates["species"] = sci_name
                elif rank_ncbi == "subspecies" and not candidates["subspecies"]: candidates["subspecies"] = sci_name
                elif rank_ncbi == "genus" and not candidates["genus"]: candidates["genus"] = sci_name
                elif rank_ncbi == "family" and not candidates["family"]: candidates["family"] = sci_name
                elif not candidates["other_valid"]: candidates["other_valid"] = sci_name
        
        if candidates["species"]: scientific_name_suggestion = candidates["species"]
        elif candidates["subspecies"]:
            parts = candidates["subspecies"].split(); scientific_name_suggestion = " ".join(parts[:2]) if len(parts) >= 2 else candidates["subspecies"]
        elif candidates["genus"]: scientific_name_suggestion = candidates["genus"]
        elif candidates["family"]: scientific_name_suggestion = candidates["family"]
        elif candidates["other_valid"]: scientific_name_suggestion = candidates["other_valid"]
        else: print(f"  No suitable scientific name found according to prioritization.")

    except Exception as e: print(f"  Error querying NCBI Entrez: {type(e).__name__} - {e}")
    return scientific_name_suggestion

# --- Helper Function: Wikipedia Image URL (Corrected) ---
def get_wikipedia_image_url(species_name_from_user, scientific_name_from_gbif=None):
    # ... (کد کامل و صحیح این تابع از پیام قبلی اینجا قرار می‌گیرد) ...
    search_candidates=[]; clean_scientific_name=None; clean_scientific_name_for_filename=None
    if scientific_name_from_gbif: temp_clean_name=scientific_name_from_gbif.split('(')[0].strip();
    if temp_clean_name: clean_scientific_name=temp_clean_name; search_candidates.append(clean_scientific_name); clean_scientific_name_for_filename=clean_scientific_name.lower().replace(" ","_")
    user_name_for_filename=None
    if species_name_from_user: should_add_user_name=True
    if clean_scientific_name:
        if species_name_from_user.lower()==clean_scientific_name.lower(): should_add_user_name=False
    if should_add_user_name: search_candidates.append(species_name_from_user)
    user_name_for_filename=species_name_from_user.lower().replace(" ","_")
    if not search_candidates: print("[WIKI_IMG] No search terms provided."); return None
    print(f"[WIKI_IMG] Attempting image for candidates: {search_candidates}"); wikipedia.set_lang("en")
    avoid_keywords_in_filename=["map","range","distribution","locator","chart","diagram","logo","icon","disambig","sound","audio","timeline","scale","reconstruction","skeleton","skull","footprint","tracks","scat","phylogeny","cladogram","taxonomy","taxobox"]
    priority_keywords=[]
    if clean_scientific_name_for_filename: priority_keywords.append(clean_scientific_name_for_filename)
    if "_" in clean_scientific_name_for_filename: priority_keywords.append(clean_scientific_name_for_filename.split("_")[0])
    if user_name_for_filename:
        if not (clean_scientific_name_for_filename and user_name_for_filename==clean_scientific_name_for_filename):
            if user_name_for_filename not in priority_keywords: priority_keywords.append(user_name_for_filename)
        if "_" in user_name_for_filename: 
            user_genus_equivalent = user_name_for_filename.split("_")[0] # تعریف فقط اگر لازم بود
            add_user_genus = True
            if clean_scientific_name_for_filename and "_" in clean_scientific_name_for_filename:
                if user_genus_equivalent == clean_scientific_name_for_filename.split("_")[0]: add_user_genus = False
            elif clean_scientific_name_for_filename and user_genus_equivalent == clean_scientific_name_for_filename: add_user_genus = False
            if add_user_genus and user_genus_equivalent not in priority_keywords: priority_keywords.append(user_genus_equivalent)
    priority_keywords=list(filter(None,dict.fromkeys(priority_keywords))); print(f"[WIKI_IMG] Priority keywords for image filename: {priority_keywords}")
    processed_search_terms=set()
    while search_candidates:
        term_to_search=search_candidates.pop(0)
        if not term_to_search or term_to_search in processed_search_terms: continue
        processed_search_terms.add(term_to_search); print(f"[WIKI_IMG] Trying Wikipedia search for term: '{term_to_search}'")
        try: 
            wiki_page = None; page_title_found = None
            try: wiki_page = wikipedia.page(term_to_search, auto_suggest=True, redirect=True); page_title_found = wiki_page.title; print(f"[WIKI_IMG] Found page directly: '{page_title_found}' for '{term_to_search}'")
            except wikipedia.exceptions.PageError:
                print(f"[WIKI_IMG] Page not found directly for '{term_to_search}'. Trying wikipedia.search()."); search_results = wikipedia.search(term_to_search, results=1)
                if search_results:
                    page_title_to_get = search_results[0]
                    if page_title_to_get in processed_search_terms: print(f"[WIKI_IMG] Page '{page_title_to_get}' already processed, skipping."); continue 
                    print(f"[WIKI_IMG] Found page via search: '{page_title_to_get}' for '{term_to_search}'"); wiki_page = wikipedia.page(page_title_to_get, auto_suggest=False, redirect=True); page_title_found = wiki_page.title
                else: print(f"[WIKI_IMG] No search results in Wikipedia for: '{term_to_search}'"); continue 
            except wikipedia.exceptions.DisambiguationError as e:
                print(f"[WIKI_IMG] Disambiguation for '{term_to_search}'. Options: {e.options[:3]}")
                if e.options: new_candidate = e.options[0]
                if new_candidate not in processed_search_terms and new_candidate not in search_candidates: search_candidates.append(new_candidate); print(f"[WIKI_IMG] Added disambiguation option '{new_candidate}' to search candidates.")
                continue 
            if wiki_page and wiki_page.images:
                print(f"[WIKI_IMG] Page '{page_title_found}' images (up to 5): {wiki_page.images[:5]}"); candidate_images_with_scores=[]
                for img_url in wiki_page.images: 
                    img_url_lower=img_url.lower()
                    if not any(ext in img_url_lower for ext in['.png','.jpg','.jpeg']): continue
                    if any(keyword in img_url_lower for keyword in avoid_keywords_in_filename): continue
                    score=0
                    for pk_word in priority_keywords:
                        if pk_word in img_url_lower: score+=5; filename_part=img_url_lower.split('/')[-1]
                        if filename_part.startswith(pk_word): score+=3
                        if clean_scientific_name_for_filename and pk_word==clean_scientific_name_for_filename and pk_word in filename_part : score+=5
                    if img_url_lower.endswith('.svg'): score-=1; 
                    candidate_images_with_scores.append({'url':img_url,'score':score})
                if not candidate_images_with_scores: print(f"[WIKI_IMG] No images passed filter for '{page_title_found}'"); continue
                sorted_images=sorted(candidate_images_with_scores,key=lambda x:x['score'],reverse=True); print(f"[WIKI_IMG] Sorted suitable images (top 3 with scores): {[{'url':i['url'][-50:],'score':i['score']} for i in sorted_images[:3]]}")
                if sorted_images and sorted_images[0]['score'] > 0: 
                    best_image_url=sorted_images[0]['url']
                    if best_image_url.startswith("//"): best_image_url="https:"+best_image_url
                    print(f"[WIKI_IMG] Best image found: {best_image_url} with score {sorted_images[0]['score']}")
                    return best_image_url 
                else: print(f"[WIKI_IMG] No image found with positive score for '{page_title_found}'.")
            elif wiki_page: print(f"[WIKI_IMG] No images listed on Wikipedia page: '{page_title_found}'")
        except wikipedia.exceptions.PageError: print(f"[WIKI_IMG] Wikipedia PageError (likely after search/disambiguation) for term: '{term_to_search}'"); continue 
        except Exception as e: print(f"[WIKI_IMG] Generic error during processing for '{term_to_search}': {str(e)}"); traceback.print_exc(); continue 
    print(f"[WIKI_IMG] No suitable Wikipedia image URL after all attempts."); return None

# --- Endpoint: Suggest Name (فعال شده) ---
@app.route('/api/suggest_name', methods=['GET'])
def suggest_name_endpoint():
    query = request.args.get('query') # انتظار نام انگلیسی
    common_headers = {'Access-Control-Allow-Origin': '*'}
    if not query:
        return jsonify({"error": "Query parameter (English common name) is required"}), 400, common_headers
    
    # فراخوانی تابع NCBI برای گرفتن پیشنهاد نام علمی
    suggestion = get_best_ncbi_suggestion_flexible(query)

    response_data = {"query": query, "scientific_name_suggestion": suggestion}
    if suggestion:
        print(f"[SUGGEST_API] Suggestion for '{query}': {suggestion}")
        return jsonify(response_data), 200, common_headers
    else:
        print(f"[SUGGEST_API] No suggestion found for '{query}'")
        response_data["message"] = f"Could not find a scientific name suggestion for '{query}' via NCBI."
        return jsonify(response_data), 404, common_headers


# --- Endpoint: Main Handler (Taxonomy & Image - بدون تغییر) ---
@app.route('/', defaults={'path': ''}, methods=['GET', 'POST', 'OPTIONS'])
@app.route('/<path:path>', methods=['GET', 'POST', 'OPTIONS'])
def main_handler(path=None):
    common_headers = {'Access-Control-Allow-Origin': '*'}
    if request.method == 'OPTIONS':
        cors_headers = {**common_headers, 'Access-Control-Allow-Methods': 'GET, POST, OPTIONS', 'Access-Control-Allow-Headers': 'Content-Type, Authorization', 'Access-Control-Max-Age': '3600'}
        return ('', 204, cors_headers)

    scientific_name_query = ""
    if request.method == 'GET': scientific_name_query = request.args.get('name')
    elif request.method == 'POST':
        try: data_post = request.get_json(); scientific_name_query = data_post.get('name') if data_post else None
        except Exception as e_post: print(f"Error parsing JSON body: {str(e_post)}"); return jsonify({"error": "Invalid JSON"}), 400, common_headers
        if not scientific_name_query: return jsonify({"error": "JSON body must contain 'name' field."}), 400, common_headers
    
    if not scientific_name_query: return jsonify({"error": "Parameter 'name' (exact scientific name) is required."}), 400, common_headers

    print(f"[MAIN_HANDLER] Processing request for scientific name: '{scientific_name_query}'")
    params_gbif = {"name": scientific_name_query, "verbose": "true"}
    classification_data = {"searchedName": scientific_name_query} 
    gbif_error_message = None
    gbif_scientific_name_confirmed = None 

    try:
        api_response_gbif = requests.get(GBIF_API_URL_MATCH, params=params_gbif, timeout=10)
        api_response_gbif.raise_for_status()
        data_gbif = api_response_gbif.json()
        if data_gbif.get("matchType") == "EXACT" and data_gbif.get("confidence", 0) > 90 :
            gbif_scientific_name_confirmed = data_gbif.get("scientificName")
            classification_data.update({"scientificName": gbif_scientific_name_confirmed, "kingdom": data_gbif.get("kingdom"), "phylum": data_gbif.get("phylum"), "class": data_gbif.get("class"), "order": data_gbif.get("order"), "family": data_gbif.get("family"), "genus": data_gbif.get("genus"), "species": data_gbif.get("species") if data_gbif.get("speciesKey") and data_gbif.get("species") else None, "usageKey": data_gbif.get("usageKey"), "confidence": data_gbif.get("confidence"), "matchType": data_gbif.get("matchType"), "status": data_gbif.get("status"), "rank": data_gbif.get("rank")})   
        else:
             gbif_error_message = f"GBIF could not find an exact match for '{scientific_name_query}'. MatchType: {data_gbif.get('matchType')}, Confidence: {data_gbif.get('confidence')}"
             print(f"[GBIF_WARN] {gbif_error_message}")
             classification_data.update({k: v for k, v in data_gbif.items() if k in ['scientificName', 'kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'usageKey', 'confidence', 'matchType', 'status', 'rank']})
    except requests.exceptions.Timeout: gbif_error_message = "GBIF Timeout"; print(f"[GBIF_ERR] Timeout for: {scientific_name_query}")
    except requests.exceptions.HTTPError as http_err_gbif: gbif_error_message = f"GBIF Server Error: {http_err_gbif}"; print(f"[GBIF_ERR] HTTPError: {gbif_error_message}")
    except requests.exceptions.RequestException as e_gbif_req: gbif_error_message = f"GBIF Network Error: {str(e_gbif_req)}"; print(f"[GBIF_ERR] RequestException: {str(e_gbif_req)}")
    except Exception as e_gbif_generic: gbif_error_message = "GBIF Internal Error."; print(f"[GBIF_ERR] Generic Error: {str(e_gbif_generic)}"); traceback.print_exc()

    name_for_wiki_image = gbif_scientific_name_confirmed if gbif_scientific_name_confirmed else scientific_name_query
    wiki_image_url = get_wikipedia_image_url(scientific_name_query, name_for_wiki_image) 
    if wiki_image_url: classification_data["imageUrl"] = wiki_image_url
    if gbif_error_message and not gbif_scientific_name_confirmed: classification_data["gbifLookupMessage"] = gbif_error_message
    final_data = {k: v for k, v in classification_data.items() if v is not None}
    return jsonify(final_data), 200, common_headers

# --- اجرای محلی (برای تست) ---
# if __name__ == "__main__":
#     app.run(debug=True, port=5001)
