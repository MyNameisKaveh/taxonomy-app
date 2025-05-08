# api/utils.py
from Bio import Entrez
import time
import traceback # اگر از traceback در تابع استفاده کردی

# !!! ایمیل خود را اینجا تنظیم کنید !!!
Entrez.email = "andolini1889@gmail.com" 

def get_best_ncbi_suggestion_flexible(common_name, max_ids_to_check=5):
    # ... (کد کامل تابع که در Colab تست شد) ...
    print(f"\n--- Processing '{common_name}' with NCBI Entrez (Flexible Ranks) ---")
    scientific_name_suggestion = None
    try:
        search_term = f"{common_name}[Common Name] OR {common_name}[Organism]"
        # ... (بقیه کد تابع) ...
        handle = Entrez.esearch(db="taxonomy", term=search_term, retmax=max_ids_to_check, sort="relevance")
        record = Entrez.read(handle)
        handle.close()
        id_list = record["IdList"]

        if not id_list:
            time.sleep(0.34)
            handle = Entrez.esearch(db="taxonomy", term=common_name, retmax=max_ids_to_check, sort="relevance")
            record = Entrez.read(handle)
            handle.close()
            id_list = record["IdList"]
            if not id_list:
                print(f"  No TaxIDs found for '{common_name}' even after retry.")
                return None
        
        print(f"  Found {len(id_list)} potential TaxIDs: {id_list}")
        ids_to_fetch = id_list[:max_ids_to_check] 
        if not ids_to_fetch: return None

        wait_time = 0.2 + (0.34 * len(ids_to_fetch))
        time.sleep(wait_time) 
        
        summary_handle = Entrez.esummary(db="taxonomy", id=",".join(ids_to_fetch), retmode="xml")
        summary_records = Entrez.read(summary_handle)
        summary_handle.close()
        
        candidates = { "species": None, "subspecies": None, "genus": None, "family": None, "other_valid": None }
        
        for summary_record in summary_records:
            sci_name = summary_record.get("ScientificName")
            rank_ncbi = summary_record.get("Rank", "N/A").lower()
            tax_id_processed = summary_record.get("Id") 

            print(f"    - Candidate TaxID {tax_id_processed}: '{sci_name}', Rank: {rank_ncbi}")

            if sci_name:
                if rank_ncbi == "species":
                    if not candidates["species"]:
                        print(f"      -> Found Species Match: '{sci_name}'")
                        candidates["species"] = sci_name
                elif rank_ncbi == "subspecies":
                    if not candidates["subspecies"]:
                        print(f"      -> Found Subspecies Match: '{sci_name}'")
                        candidates["subspecies"] = sci_name
                elif rank_ncbi == "genus":
                     if not candidates["genus"]:
                         print(f"      -> Found Genus Match: '{sci_name}'")
                         candidates["genus"] = sci_name
                elif rank_ncbi == "family":
                     if not candidates["family"]:
                         print(f"      -> Found Family Match: '{sci_name}'")
                         candidates["family"] = sci_name
                elif not candidates["other_valid"]:
                     print(f"      -> Found Other Valid Scientific Name: '{sci_name}' (Rank: {rank_ncbi})")
                     candidates["other_valid"] = sci_name
        
        if candidates["species"]:
            scientific_name_suggestion = candidates["species"]
        elif candidates["subspecies"]:
            parts = candidates["subspecies"].split()
            if len(parts) >= 2:
                species_name_from_subspecies = " ".join(parts[:2])
                print(f"  -> Using Species name from Subspecies: '{species_name_from_subspecies}'")
                scientific_name_suggestion = species_name_from_subspecies
            else: 
                scientific_name_suggestion = candidates["subspecies"] 
        elif candidates["genus"]:
            print(f"  -> Using Genus name as best match: '{candidates['genus']}'")
            scientific_name_suggestion = candidates["genus"]
        elif candidates["family"]:
            print(f"  -> Using Family name as best match: '{candidates['family']}'")
            scientific_name_suggestion = candidates["family"]
        elif candidates["other_valid"]:
             print(f"  -> Using first other valid name as best match: '{candidates['other_valid']}'")
             scientific_name_suggestion = candidates["other_valid"]
        else:
            print(f"  No suitable scientific name found according to prioritization.")

    except Exception as e:
        print(f"  Error querying NCBI Entrez: {type(e).__name__} - {e}")
        # traceback.print_exc() # برای دیباگ بیشتر فعال شود
        
    return scientific_name_suggestion

# نکته: مطمئن شوید تابع get_wikipedia_image_url هم در utils.py قرار دارد یا در handler.py باقی می‌ماند و import می‌شود.
# برای سادگی، فرض می‌کنیم تابع تصویر در handler.py باقی می‌ماند.
