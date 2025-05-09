# تابع به‌روز شده برای Colab

from Bio import Entrez
import time
import json

Entrez.email = "YOUR_ACTUAL_EMAIL@example.com" # !!! ایمیل خودتون رو اینجا بذارید !!!

def get_best_ncbi_suggestion_flexible(common_name, max_ids_to_check=5): # افزایش تعداد برای شانس بیشتر
    """
    با استفاده از NCBI Entrez، بهترین نام علمی ممکن را برای یک نام رایج پیدا می‌کند،
    با اولویت‌بندی رتبه‌های مختلف (گونه > زیرگونه > جنس > خانواده > سایر).
    """
    print(f"\n--- Processing '{common_name}' with NCBI Entrez (Flexible Ranks) ---")
    scientific_name_suggestion = None
    try:
        search_term = f"{common_name}[Common Name] OR {common_name}[Organism]"
        # print(f"  Searching NCBI Taxonomy with term: '{search_term}'")
        handle = Entrez.esearch(db="taxonomy", term=search_term, retmax=max_ids_to_check, sort="relevance")
        record = Entrez.read(handle)
        handle.close()
        id_list = record["IdList"]

        if not id_list:
            # print(f"  No TaxIDs found for '{common_name}' with field filter. Retrying without filter...")
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

        # افزایش تاخیر متناسب با تعداد ID ها
        wait_time = 0.2 + (0.34 * len(ids_to_fetch)) # حداقل 0.5 ثانیه
        # print(f"  Waiting for {wait_time:.2f} seconds before esummary...")
        time.sleep(wait_time) 
        
        summary_handle = Entrez.esummary(db="taxonomy", id=",".join(ids_to_fetch), retmode="xml")
        summary_records = Entrez.read(summary_handle)
        summary_handle.close()
        
        # ذخیره بهترین کاندید برای هر سطح اولویت
        candidates = { 
            "species": None, 
            "subspecies": None, # خود نام زیرگونه رو نگه میداریم
            "genus": None, 
            "family": None,
            "other_valid": None # اولین نام علمی معتبر با هر رتبه دیگر
        }
        
        for summary_record in summary_records:
            sci_name = summary_record.get("ScientificName")
            rank_ncbi = summary_record.get("Rank", "N/A").lower()
            tax_id_processed = summary_record.get("Id") 

            print(f"    - Candidate TaxID {tax_id_processed}: '{sci_name}', Rank: {rank_ncbi}")

            if sci_name:
                if rank_ncbi == "species":
                    if not candidates["species"]: # اولین گونه پیدا شده
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
                # اگر هیچکدام از موارد بالا نبود و اولین نام علمی معتبر بود
                elif not candidates["other_valid"]:
                     print(f"      -> Found Other Valid Scientific Name: '{sci_name}' (Rank: {rank_ncbi})")
                     candidates["other_valid"] = sci_name
        
        # تصمیم‌گیری نهایی بر اساس اولویت
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
        
    return scientific_name_suggestion

# تست‌ها با نام‌های متنوع‌تر
test_common_names_flexible = ["Rose", "Lion", "Human", "Cat", "Apple tree", "Escherichia coli", "Grizzly Bear", "Blue Whale", "Ladybug", "Oak", "Felidae", "Bacteria"]
final_suggestions_flexible = {}
for name in test_common_names_flexible:
    suggestion = get_best_ncbi_suggestion_flexible(name)
    final_suggestions_flexible[name] = suggestion
    print(f"Final Flexible NCBI Suggestion for '{name}': {suggestion}\n")

# نمایش خلاصه نتایج
print("\n--- Summary of Flexible NCBI Suggestions ---")
for name, suggestion in final_suggestions_flexible.items():
    print(f"'{name}' => {suggestion}")
