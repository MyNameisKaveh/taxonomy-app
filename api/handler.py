from flask import Flask, jsonify, request
import requests
import traceback
import wikipedia
import time # برای NCBI

# --- NCBI Entrez Imports and Setup ---
from Bio import Entrez
Entrez.email = "YOUR_ACTUAL_EMAIL@example.com" # !!! ایمیل معتبر خودتان را اینجا وارد کنید !!!
# --- END NCBI ---

app = Flask(__name__)
GBIF_API_URL_MATCH = "https://api.gbif.org/v1/species/match"

# --- START: NCBI Suggestion Function (کد تست شما با کمی تغییرات جزئی برای ادغام) ---
def get_best_ncbi_suggestion_flexible(common_name, max_ids_to_check=5):
    app.logger.info(f"[NCBI_LOG] Processing '{common_name}' with NCBI Entrez (Flexible Ranks)")
    scientific_name_suggestion = None
    try:
        search_term = f"{common_name}[Common Name] OR {common_name}[Organism]"
        handle = Entrez.esearch(db="taxonomy", term=search_term, retmax=max_ids_to_check, sort="relevance")
        record = Entrez.read(handle)
        handle.close()
        id_list = record["IdList"]

        if not id_list:
            app.logger.info(f"[NCBI_LOG] No TaxIDs found for '{common_name}' with field filter. Retrying without filter...")
            time.sleep(0.34) # رعایت نرخ درخواست NCBI
            handle = Entrez.esearch(db="taxonomy", term=common_name, retmax=max_ids_to_check, sort="relevance")
            record = Entrez.read(handle)
            handle.close()
            id_list = record["IdList"]
            if not id_list:
                app.logger.info(f"[NCBI_LOG] No TaxIDs found for '{common_name}' even after retry.")
                return None
        
        app.logger.info(f"[NCBI_LOG] Found {len(id_list)} potential TaxIDs: {id_list}")
        ids_to_fetch = id_list[:max_ids_to_check] 
        if not ids_to_fetch: return None

        wait_time = 0.2 + (0.34 * len(ids_to_fetch))
        time.sleep(wait_time) 
        
        summary_handle = Entrez.esummary(db="taxonomy", id=",".join(ids_to_fetch), retmode="xml")
        summary_records = Entrez.read(summary_handle)
        summary_handle.close()
        
        candidates = { 
            "species": None, "subspecies": None, "genus": None, 
            "family": None, "other_valid": None
        }
        
        for summary_record in summary_records:
            sci_name = summary_record.get("ScientificName")
            rank_ncbi = summary_record.get("Rank", "N/A").lower()
            tax_id_processed = summary_record.get("Id") 
            app.logger.debug(f"[NCBI_LOG] Candidate TaxID {tax_id_processed}: '{sci_name}', Rank: {rank_ncbi}")

            if sci_name:
                if rank_ncbi == "species" and not candidates["species"]:
                     candidates["species"] = sci_name
                elif rank_ncbi == "subspecies" and not candidates["subspecies"]:
                     candidates["subspecies"] = sci_name
                elif rank_ncbi == "genus" and not candidates["genus"]:
                     candidates["genus"] = sci_name
                elif rank_ncbi == "family" and not candidates["family"]:
                     candidates["family"] = sci_name
                elif not candidates["other_valid"]: # اولین نام علمی معتبر دیگر
                     candidates["other_valid"] = sci_name
        
        if candidates["species"]:
            scientific_name_suggestion = candidates["species"]
        elif candidates["subspecies"]:
            parts = candidates["subspecies"].split()
            scientific_name_suggestion = " ".join(parts[:2]) if len(parts) >= 2 else candidates["subspecies"]
        elif candidates["genus"]:
            scientific_name_suggestion = candidates["genus"]
        elif candidates["family"]:
            scientific_name_suggestion = candidates["family"]
        elif candidates["other_valid"]:
             scientific_name_suggestion = candidates["other_valid"]
        else:
            app.logger.info(f"[NCBI_LOG] No suitable scientific name found according to prioritization for '{common_name}'.")

    except Exception as e:
        app.logger.error(f"[NCBI_ERR] Error querying NCBI Entrez for '{common_name}': {type(e).__name__} - {e}")
        traceback.print_exc()
        
    app.logger.info(f"[NCBI_LOG] Final NCBI Suggestion for '{common_name}': {scientific_name_suggestion}")
    return scientific_name_suggestion
# --- END: NCBI Suggestion Function ---


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
        app.logger.warning("[WIKI_IMG] No search terms provided.")
        return None

    app.logger.info(f"[WIKI_IMG] Attempting image for candidates: {search_candidates}")
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
    app.logger.info(f"[WIKI_IMG] Priority keywords for image filename: {priority_keywords}")

    processed_search_terms = set() 

    while search_candidates:
        term_to_search = search_candidates.pop(0) 
        if not term_to_search or term_to_search in processed_search_terms:
            continue
        processed_search_terms.add(term_to_search)
        
        app.logger.info(f"[WIKI_IMG] Trying Wikipedia search for term: '{term_to_search}'")
        try:
            wiki_page = None
            try:
                wiki_page = wikipedia.page(term_to_search, auto_suggest=True, redirect=True)
                app.logger.info(f"[WIKI_IMG] Found page directly: '{wiki_page.title}' for '{term_to_search}'")
            except wikipedia.exceptions.PageError:
                app.logger.info(f"[WIKI_IMG] Page not found directly for '{term_to_search}'. Trying wikipedia.search().")
                search_results = wikipedia.search(term_to_search, results=1)
                if not search_results:
                    app.logger.info(f"[WIKI_IMG] No search results in Wikipedia for: '{term_to_search}'")
                    continue
                page_title_to_get = search_results[0]
                if page_title_to_get in processed_search_terms: 
                    app.logger.info(f"[WIKI_IMG] Page '{page_title_to_get}' already processed, skipping.")
                    continue
                app.logger.info(f"[WIKI_IMG] Found page via search: '{page_title_to_get}' for '{term_to_search}'")
                wiki_page = wikipedia.page(page_title_to_get, auto_suggest=False, redirect=True)
            
            if wiki_page and wiki_page.images:
                app.logger.debug(f"[WIKI_IMG] Page '{wiki_page.title}' images (up to 5): {wiki_page.images[:5]}")
                
                candidate_images_with_scores = []
                for img_url in wiki_page.images:
                    img_url_lower = img_url.lower()
                    if not any(ext in img_url_lower for ext in ['.png', '.jpg', '.jpeg']): 
                        continue
                    if any(keyword in img_url_lower for keyword in avoid_keywords_in_filename):
                        continue
                    
                    score = 0
                    for pk_word in priority_keywords:
                        if pk_word in img_url_lower:
                            score += 5 
                            filename_part = img_url_lower.split('/')[-1]
                            if filename_part.startswith(pk_word): 
                                score += 3
                            if pk_word == clean_scientific_name_for_filename and pk_word in filename_part : 
                                score +=5
                    if img_url_lower.endswith('.svg'): score -= 1
                    candidate_images_with_scores.append({'url': img_url, 'score': score})
                
                if not candidate_images_with_scores:
                    app.logger.info(f"[WIKI_IMG] No images passed filter for '{wiki_page.title}'")
                    continue

                sorted_images = sorted(candidate_images_with_scores, key=lambda x: x['score'], reverse=True)
                app.logger.debug(f"[WIKI_IMG] Sorted suitable images (top 3 with scores): {[{'url': i['url'][-50:], 'score': i['score']} for i in sorted_images[:3]]}")

                if sorted_images and sorted_images[0]['score'] > 0: 
                    best_image_url = sorted_images[0]['url']
                    if best_image_url.startswith("//"): best_image_url = "https:" + best_image_url
                    app.logger.info(f"[WIKI_IMG] Best image found: {best_image_url} with score {sorted_images[0]['score']}")
                    return best_image_url
                else: 
                    app.logger.info(f"[WIKI_IMG] No image found with positive score for '{wiki_page.title}'.")

            elif wiki_page:
                app.logger.info(f"[WIKI_IMG] No images listed on Wikipedia page: '{wiki_page.title}'")
            
        except wikipedia.exceptions.DisambiguationError as e:
            app.logger.warning(f"[WIKI_IMG] Disambiguation for '{term_to_search}'. Options: {e.options[:3]}")
            if e.options:
                new_candidate = e.options[0]
                if new_candidate not in processed_search_terms and new_candidate not in search_candidates:
                    search_candidates.append(new_candidate)
                    app.logger.info(f"[WIKI_IMG] Added disambiguation option '{new_candidate}' to search candidates.")
            continue 
        except wikipedia.exceptions.PageError:
             app.logger.warning(f"[WIKI_IMG] Wikipedia PageError (likely after search or disambiguation) for term: '{term_to_search}'")
             continue
        except Exception as e:
            app.logger.error(f"[WIKI_IMG] Generic error for '{term_to_search}': {str(e)}")
            traceback.print_exc() 
            continue
            
    app.logger.warning(f"[WIKI_IMG] No suitable Wikipedia image URL after all attempts for initial candidates.")
    return None


def _fetch_gbif_data(name_to_search, search_source="initial"):
    """Helper function to query GBIF and parse response."""
    params_gbif = {"name": name_to_search, "verbose": "true"}
    data = {}
    error_message = None
    scientific_name_found = None
    
    app.logger.info(f"[GBIF_API] Querying GBIF for '{name_to_search}' (source: {search_source})")
    try:
        api_response_gbif = requests.get(GBIF_API_URL_MATCH, params=params_gbif, timeout=10)
        api_response_gbif.raise_for_status()
        gbif_json = api_response_gbif.json()
        
        if not gbif_json or gbif_json.get("matchType") == "NONE" or gbif_json.get("confidence", 0) < 30:
            error_message = f"موجودی با نام '{name_to_search}' در GBIF پیدا نشد یا نتیجه با اطمینان کافی نبود. (منبع جستجو: {search_source})"
            data = {"searchedName": name_to_search, "matchType": gbif_json.get("matchType", "NONE"), "confidence": gbif_json.get("confidence")}
            app.logger.warning(f"[GBIF_API] Low confidence/no match for '{name_to_search}'. MatchType: {data.get('matchType')}, Confidence: {data.get('confidence')}")
        else:
            scientific_name_found = gbif_json.get("scientificName")
            data = {
                "scientificName": scientific_name_found,
                "kingdom": gbif_json.get("kingdom"), "phylum": gbif_json.get("phylum"), "class": gbif_json.get("class"),
                "order": gbif_json.get("order"), "family": gbif_json.get("family"), "genus": gbif_json.get("genus"),
                "species": gbif_json.get("species") if gbif_json.get("speciesKey") and gbif_json.get("species") else None,
                "usageKey": gbif_json.get("usageKey"), "confidence": gbif_json.get("confidence"),
                "matchType": gbif_json.get("matchType"), "status": gbif_json.get("status"), "rank": gbif_json.get("rank")
            }
            app.logger.info(f"[GBIF_API] Success for '{name_to_search}'. ScientificName: {scientific_name_found}, Confidence: {data['confidence']}")
            
    except requests.exceptions.Timeout:
        error_message = f"خطا: زمان پاسخگویی از سرور GBIF برای '{name_to_search}' بیش از حد طول کشید."
        app.logger.error(f"[GBIF_ERR] Timeout for: {name_to_search}")
    except requests.exceptions.HTTPError as http_err_gbif:
        error_message = f"خطا از سرور GBIF برای '{name_to_search}': {http_err_gbif}"
        try:
            gbif_error_details = api_response_gbif.json() # api_response_gbif might not be defined if error is early
            error_message += f" - پیام GBIF: {gbif_error_details.get('message', api_response_gbif.text[:100])}"
        except: pass
        app.logger.error(f"[GBIF_ERR] HTTPError for '{name_to_search}': {error_message}")
    except requests.exceptions.RequestException as e_gbif_req:
        error_message = f"خطا در ارتباط با سرور GBIF برای '{name_to_search}': {str(e_gbif_req)}"
        app.logger.error(f"[GBIF_ERR] RequestException for '{name_to_search}': {str(e_gbif_req)}")
    except Exception as e_gbif_generic:
        error_message = f"یک خطای پیش‌بینی نشده داخلی در ارتباط با GBIF برای '{name_to_search}' رخ داده است."
        app.logger.error(f"[GBIF_ERR] Generic Error for '{name_to_search}': {str(e_gbif_generic)}")
        traceback.print_exc()
        
    return data, scientific_name_found, error_message


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
            app.logger.error(f"Error parsing JSON body: {str(e_post)}\n{traceback.format_exc()}")
            return jsonify({"error": f"خطا در پردازش درخواست: {str(e_post)}"}), 400, common_headers

    if not species_name_query:
        return jsonify({"error": "پارامتر 'name' (در آدرس یا بدنه JSON) مورد نیاز است."}), 400, common_headers

    app.logger.info(f"--- New request for '{species_name_query}' ---")

    # --- Initial GBIF Search ---
    classification_data, gbif_scientific_name, gbif_error_message = _fetch_gbif_data(species_name_query, search_source="initial user query")
    classification_data["searchedName"] = species_name_query # Always include original search term
    
    ncbi_suggestion_used = False
    ncbi_suggested_name_for_gbif = None

    # --- If initial GBIF failed or low confidence, try NCBI suggestion ---
    if gbif_error_message or not gbif_scientific_name: # gbif_error_message implies no scientific name
        app.logger.info(f"[GBIF_PRIMARY_FAIL] Initial GBIF search for '{species_name_query}' was not successful. Error: {gbif_error_message}. Trying NCBI suggestion.")
        
        original_gbif_error = gbif_error_message # Store original error
        original_classification_data = classification_data.copy() if classification_data else {} # Store original (failed) data

        ncbi_suggested_name_for_gbif = get_best_ncbi_suggestion_flexible(species_name_query)
        
        if ncbi_suggested_name_for_gbif:
            classification_data["ncbiSuggestedName"] = ncbi_suggested_name_for_gbif
            ncbi_suggestion_used = True
            app.logger.info(f"[NCBI_RETRY] NCBI suggested '{ncbi_suggested_name_for_gbif}'. Re-querying GBIF.")
            
            # --- GBIF Search with NCBI suggestion ---
            classification_data_ncbi, gbif_scientific_name_ncbi, gbif_error_message_ncbi = _fetch_gbif_data(ncbi_suggested_name_for_gbif, search_source="NCBI suggestion")
            
            if gbif_scientific_name_ncbi: # Success with NCBI suggestion
                app.logger.info(f"[GBIF_NCBI_SUCCESS] GBIF found match using NCBI suggestion: '{gbif_scientific_name_ncbi}'")
                classification_data = classification_data_ncbi # Use new data
                classification_data["searchedName"] = species_name_query # Ensure original search term is present
                classification_data["ncbiSuggestedName"] = ncbi_suggested_name_for_gbif # Add NCBI name
                gbif_scientific_name = gbif_scientific_name_ncbi # Update main scientific name
                gbif_error_message = None # Clear previous error
            else: # GBIF still failed even with NCBI suggestion
                app.logger.warning(f"[GBIF_NCBI_FAIL] GBIF search failed even with NCBI suggestion '{ncbi_suggested_name_for_gbif}'. Error: {gbif_error_message_ncbi}")
                # Revert to original error or use NCBI-related error if more specific
                gbif_error_message = f"جستجوی اولیه ناموفق بود. پیشنهاد NCBI ({ncbi_suggested_name_for_gbif}) نیز نتیجه‌ای در GBIF نداشت. خطای GBIF برای پیشنهاد NCBI: {gbif_error_message_ncbi or 'نامشخص'}"
                classification_data = original_classification_data # Use original failed data
                classification_data["searchedName"] = species_name_query
                classification_data["ncbiSuggestedName"] = ncbi_suggested_name_for_gbif
                classification_data["ncbiSuggestionUsedButGIBFFailed"] = True
                gbif_scientific_name = None # Ensure it's None
        else: # NCBI found no suggestion
            app.logger.info(f"[NCBI_NO_SUGGESTION] NCBI found no suggestion for '{species_name_query}'. Using original GBIF error.")
            gbif_error_message = original_gbif_error # Use the original GBIF error
            classification_data = original_classification_data
            classification_data["searchedName"] = species_name_query
            gbif_scientific_name = None
            
    # --- Wikipedia Image Search ---
    # Uses gbif_scientific_name which might be from the NCBI-assisted GBIF search
    wiki_image_url = get_wikipedia_image_url(species_name_query, gbif_scientific_name) 
    if wiki_image_url:
        classification_data["imageUrl"] = wiki_image_url
    
    # --- Prepare Final Response ---
    final_data = {k: v for k, v in classification_data.items() if v is not None}
    if ncbi_suggestion_used:
        final_data["ncbiSuggestionUsed"] = True
        if ncbi_suggested_name_for_gbif: # Should always be true if ncbi_suggestion_used is true
             final_data["ncbiOriginalSuggestion"] = ncbi_suggested_name_for_gbif

    if gbif_error_message and not final_data.get("scientificName"): 
        # If there's an error AND no scientific name was ultimately found
        response_payload = {"message": gbif_error_message, **final_data} # final_data might have 'searchedName', 'ncbiSuggestedName', 'imageUrl' etc.
        if not final_data.get("imageUrl"): # If also no image, it's a more definite "not found"
            response_payload["searchedNameOnly"] = species_name_query # ensure searchedName is there if final_data was empty
        return jsonify(response_payload), 404, common_headers
    
    # If we have a scientific name, or even just an image, consider it a partial/full success
    return jsonify(final_data), 200, common_headers

if __name__ == '__main__':
    # برای تست محلی بهتر است از logger فلاسک استفاده شود
    # app.debug = True # برای نمایش لاگ‌های دیباگ بیشتر در کنسول
    # import logging
    # logging.basicConfig(level=logging.DEBUG) # این هم به نمایش لاگ‌ها کمک می‌کند
    # handler = logging.StreamHandler()
    # handler.setLevel(logging.DEBUG)
    # app.logger.addHandler(handler)
    app.run(host='0.0.0.0', port=5000) # در محیط پروداکشن از Gunicorn یا مشابه استفاده کنید
