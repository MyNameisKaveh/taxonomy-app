// URLs for API endpoints
const API_BASE_URL_MAIN_SEARCH = "https://taxonomy-app-ebon.vercel.app/api/handler"; 
const API_BASE_URL_SUGGESTION = "https://taxonomy-app-ebon.vercel.app/api/suggest_scientific_name"; 
const GBIF_SUGGEST_API_URL = "https://api.gbif.org/v1/species/suggest";

// Elements for Main Search
const speciesNameInput = document.getElementById('speciesNameInput');
const searchButton = document.getElementById('searchButton');
const resultsContainer = document.getElementById('resultsContainer');
const loadingIndicator = document.getElementById('loadingIndicator');
const errorContainer = document.getElementById('errorContainer');
const speciesImageContainer = document.getElementById('speciesImageContainer');
const autocompleteSuggestionsContainer = document.getElementById('autocompleteSuggestions'); // New element

// Elements for Suggestion Feature
const commonNameInput = document.getElementById('commonNameInput');
const suggestNameButton = document.getElementById('suggestNameButton');
const suggestionResultContainer = document.getElementById('suggestionResultContainer');
const suggestionLoadingIndicator = document.getElementById('suggestionLoadingIndicator');

// --- Debounce Function ---
function debounce(func, delay) {
    let timeoutId;
    return function(...args) {
        clearTimeout(timeoutId);
        timeoutId = setTimeout(() => {
            func.apply(this, args);
        }, delay);
    };
}

// --- Event Listeners for Main Search ---
if (searchButton) {
    searchButton.addEventListener('click', performMainSearch);
}
if (speciesNameInput) {
    speciesNameInput.addEventListener('keypress', function(event) {
        if (event.key === 'Enter') {
            // If suggestions are visible, the first one might be selected by Enter.
            // Or, we might want to explicitly hide suggestions and then search.
            // For now, just perform search. If a suggestion was clicked, it would have filled the input.
            if (autocompleteSuggestionsContainer.style.display === 'block' && autocompleteSuggestionsContainer.firstChild) {
                 // Potentially trigger click on first suggestion or handle differently.
                 // For now, let's hide suggestions and let the current input value be searched.
                autocompleteSuggestionsContainer.innerHTML = '';
                autocompleteSuggestionsContainer.style.display = 'none';
            }
            performMainSearch();
        }
    });

    // Autocomplete listener
    speciesNameInput.addEventListener('input', debounce(handleAutocompleteInput, 300));

    speciesNameInput.addEventListener('blur', function() {
        // Delay hiding to allow click on suggestion item
        setTimeout(() => {
            if (!autocompleteSuggestionsContainer.matches(':hover')) { // Don't hide if mouse is over suggestions
                autocompleteSuggestionsContainer.innerHTML = '';
                autocompleteSuggestionsContainer.style.display = 'none';
            }
        }, 150);
    });
}

// --- Event Listeners for Suggestion Feature ---
if (suggestNameButton) {
    suggestNameButton.addEventListener('click', performSuggestionSearch);
}
if (commonNameInput) {
    commonNameInput.addEventListener('keypress', function(event) {
        if (event.key === 'Enter') {
            performSuggestionSearch();
        }
    });
}

// --- Autocomplete Logic for Main Search Input ---
async function handleAutocompleteInput() {
    const trimmedQuery = speciesNameInput.value.trim();

    if (trimmedQuery.length < 3) {
        autocompleteSuggestionsContainer.innerHTML = '';
        autocompleteSuggestionsContainer.style.display = 'none';
        return;
    }

    autocompleteSuggestionsContainer.innerHTML = '<div>درحال جستجو...</div>'; // Loading state
    autocompleteSuggestionsContainer.style.display = 'block';

    try {
        const response = await fetch(`${GBIF_SUGGEST_API_URL}?q=${encodeURIComponent(trimmedQuery)}&limit=7`); // Limit to 7 suggestions
        if (!response.ok) {
            throw new Error(`GBIF API error: ${response.status}`);
        }
        const suggestions = await response.json();

        autocompleteSuggestionsContainer.innerHTML = ''; // Clear loading state

        if (suggestions && suggestions.length > 0) {
            suggestions.forEach(suggestion => {
                const itemDiv = document.createElement('div');
                // Prioritize scientificName, fallback to canonicalName or other relevant fields
                const nameToDisplay = suggestion.scientificName || suggestion.canonicalName || suggestion.vernacularName;
                let displayText = nameToDisplay;

                if (suggestion.scientificName && suggestion.vernacularName) {
                    displayText = `<em>${suggestion.scientificName}</em> (${suggestion.vernacularName})`;
                } else if (suggestion.canonicalName && suggestion.vernacularName) {
                    displayText = `<em>${suggestion.canonicalName}</em> (${suggestion.vernacularName})`;
                } else if (suggestion.scientificName) {
                     displayText = `<em>${suggestion.scientificName}</em>`;
                } else if (suggestion.canonicalName) {
                     displayText = `<em>${suggestion.canonicalName}</em>`;
                }


                if (nameToDisplay) { // Only add if we have something to display
                    itemDiv.innerHTML = displayText; // Use innerHTML for em tags
                    itemDiv.addEventListener('click', function() {
                        speciesNameInput.value = suggestion.scientificName || suggestion.canonicalName; // Use the more precise name for search
                        autocompleteSuggestionsContainer.innerHTML = '';
                        autocompleteSuggestionsContainer.style.display = 'none';
                        performMainSearch(); // Optionally trigger search automatically
                    });
                    autocompleteSuggestionsContainer.appendChild(itemDiv);
                }
            });
            if (autocompleteSuggestionsContainer.children.length > 0) {
                 autocompleteSuggestionsContainer.style.display = 'block';
            } else {
                autocompleteSuggestionsContainer.innerHTML = '<div>نتیجه‌ای یافت نشد.</div>';
                // Keep it displayed for a short while or hide immediately
                // setTimeout(() => { autocompleteSuggestionsContainer.style.display = 'none'; }, 2000);
            }
        } else {
            autocompleteSuggestionsContainer.innerHTML = '<div>نتیجه‌ای یافت نشد.</div>';
            // Keep it displayed for a short while or hide immediately
            // setTimeout(() => { autocompleteSuggestionsContainer.style.display = 'none'; }, 2000);
        }
    } catch (error) {
        console.error("Autocomplete fetch error:", error);
        autocompleteSuggestionsContainer.innerHTML = '<div>خطا در دریافت پیشنهادات.</div>';
        // Keep error displayed for a bit or hide
        // setTimeout(() => { autocompleteSuggestionsContainer.style.display = 'none'; }, 2000);
    }
    // Ensure it's hidden if no children were appended (e.g. error or no results and not keeping message displayed)
    if (autocompleteSuggestionsContainer.children.length === 0) {
        autocompleteSuggestionsContainer.style.display = 'none';
    }
}

// Global click listener to hide autocomplete suggestions
document.addEventListener('click', function(event) {
    if (speciesNameInput && autocompleteSuggestionsContainer) {
        const isClickInsideInput = speciesNameInput.contains(event.target);
        const isClickInsideSuggestions = autocompleteSuggestionsContainer.contains(event.target);
        if (!isClickInsideInput && !isClickInsideSuggestions) {
            autocompleteSuggestionsContainer.innerHTML = '';
            autocompleteSuggestionsContainer.style.display = 'none';
        }
    }
});


// --- Function to Perform Suggestion Search (NCBI) ---
async function performSuggestionSearch() {
    const commonName = commonNameInput.value.trim();
    
    suggestionResultContainer.innerHTML = '';
    suggestionResultContainer.style.display = 'none';
    suggestionLoadingIndicator.style.display = 'none';

    if (!commonName) {
        showSuggestionMessage("لطفاً نام رایج یک موجود را برای پیشنهاد وارد کنید.", "error");
        return;
    }

    suggestionLoadingIndicator.style.display = 'block'; 

    try {
        const suggestApiUrl = `${API_BASE_URL_SUGGESTION}?common_name=${encodeURIComponent(commonName)}`;
        const response = await fetch(suggestApiUrl);
        const data = await response.json();

        suggestionLoadingIndicator.style.display = 'none'; 

        if (!response.ok) {
            showSuggestionMessage(
                data.error || data.message || `مشکلی در ارتباط با سرویس پیشنهاد نام رخ داد (کد: ${response.status}). لطفاً دوباره تلاش کنید.`, 
                "error"
            ); 
        } else {
            if (data.suggested_scientific_name) {
                showSuggestionMessage(
                    `نام علمی پیشنهادی برای '${data.common_name_searched}': <br><strong>${data.suggested_scientific_name}</strong>
                     <button class="copy-to-main-search-btn" data-name="${data.suggested_scientific_name}">استفاده از این نام</button>`, 
                    "success"
                );
                const copyBtn = suggestionResultContainer.querySelector('.copy-to-main-search-btn');
                if (copyBtn) {
                    copyBtn.addEventListener('click', function() {
                        speciesNameInput.value = this.dataset.name;
                        speciesNameInput.focus();
                        showSuggestionMessage(
                            `نام علمی '${this.dataset.name}' به کادر جستجوی اصلی منتقل شد. اکنون می‌توانید دکمه "جستجو" را بزنید.`, 
                            "info"
                        ); 
                        speciesNameInput.scrollIntoView({ behavior: 'smooth', block: 'center' });
                    });
                }
            } else {
                showSuggestionMessage(
                    data.message || `متاسفانه، نام علمی دقیقی برای '${commonName}' از طریق NCBI یافت نشد. لطفاً نام دیگری را امتحان کنید یا از نام علمی شناخته شده استفاده نمایید.`, 
                    "info"
                ); 
            }
        }

    } catch (error) {
        suggestionLoadingIndicator.style.display = 'none'; 
        showSuggestionMessage(`مشکلی در ارتباط با سرویس پیشنهاد نام رخ داد. لطفاً اتصال اینترنت خود را بررسی کرده و دوباره تلاش کنید. (پیام سیستم: ${error.message})`, "error"); 
        console.error("Suggestion Fetch Error:", error);
    }
}

// --- Function to Display Messages for Suggestion Box ---
function showSuggestionMessage(message, type = "info") { 
    let cssClass = 'suggestion-message-info'; 
    if (type === 'success') {
        cssClass = 'suggested-name-success';
    } else if (type === 'error') {
        cssClass = 'suggestion-message-error';
    }
    suggestionResultContainer.innerHTML = `<p class="${cssClass}">${message}</p>`;
    suggestionResultContainer.style.display = 'block';
}


// --- Function to Perform Main Search (GBIF & Wikipedia) ---
async function performMainSearch() {
    const speciesName = speciesNameInput.value.trim();

    // Hide autocomplete when search is performed
    if (autocompleteSuggestionsContainer) {
        autocompleteSuggestionsContainer.innerHTML = '';
        autocompleteSuggestionsContainer.style.display = 'none';
    }

    resultsContainer.innerHTML = ''; 
    resultsContainer.style.display = 'none';
    errorContainer.innerHTML = '';   
    errorContainer.style.display = 'none'; 
    if (speciesImageContainer) {
        speciesImageContainer.innerHTML = ''; 
    }
    loadingIndicator.style.display = 'block'; 

    if (!speciesName) {
        showMainSearchError("لطفاً نام یک موجود را برای جستجوی اطلاعات وارد کنید."); 
        return;
    }

    try {
        const apiUrl = `${API_BASE_URL_MAIN_SEARCH}?name=${encodeURIComponent(speciesName)}`;
        const response = await fetch(apiUrl);
        const data = await response.json();

        loadingIndicator.style.display = 'none'; 

        if (!response.ok) {
            let errorMessage = data.error || data.message || `خطای ناشناخته از سرور اصلی (کد: ${response.status})`;
            if (response.status === 404 && data.message && !data.error) {
                if (data.imageUrl && speciesImageContainer) {
                    displayImageForMainSearch(data); 
                    showMainSearchInfo(data.message || `اطلاعات طبقه‌بندی برای '${speciesName}' در GBIF یافت نشد یا دقت کافی نداشت. اگر نام رایج وارد کرده‌اید، سعی کنید نام علمی آن را پیدا و جستجو کنید.`); 
                } else {
                    showMainSearchInfo(data.message || `اطلاعات طبقه‌بندی برای '${speciesName}' در GBIF یافت نشد یا دقت کافی نداشت. اگر نام رایج وارد کرده‌اید، سعی کنید نام علمی آن را پیدا و جستجو کنید.`); 
                }
            } else {
                showMainSearchError(errorMessage); 
            }
        } else {
            if (data.imageUrl && speciesImageContainer) {
                speciesImageContainer.innerHTML = '<div class="loader"></div>'; 
            } else if (speciesImageContainer) {
                 speciesImageContainer.innerHTML = ''; 
            }
            displayMainSearchResults(data);
        }

    } catch (error) {
        showMainSearchError(`مشکلی در ارتباط با سرور اصلی رخ داد. لطفاً اتصال اینترنت خود را بررسی کرده و دوباره تلاش کنید. (پیام سیستم: ${error.message})`); 
        console.error("Main Search Fetch Error:", error);
    }
}

// --- Function to Display Image for Main Search ---
function displayImageForMainSearch(data) {
    if (speciesImageContainer) {
        if (data.imageUrl) {
            const img = new Image();
            img.onload = function() {
                speciesImageContainer.innerHTML = ''; 
                speciesImageContainer.appendChild(img);
            };
            img.onerror = function() {
                speciesImageContainer.innerHTML = '<p style="text-align:center; color: #721c24;">خطا در بارگذاری تصویر. ممکن است آدرس تصویر نامعتبر باشد.</p>'; 
            };
            img.src = data.imageUrl;
            img.alt = data.scientificName || data.searchedName || 'تصویر موجود';
            img.style.maxWidth = "100%";
            img.style.maxHeight = "350px";
            img.style.borderRadius = "8px";
            img.style.display = "block"; 
            img.style.margin = "0 auto 15px auto"; 
            img.style.boxShadow = "0 4px 8px rgba(0,0,0,0.1)";
        } else {
            if (speciesImageContainer.innerHTML.includes('loader')) {
                 speciesImageContainer.innerHTML = '<p style="text-align:center; color: #7f8c8d;">تصویری برای این موجود یافت نشد.</p>';
            } else {
                 speciesImageContainer.innerHTML = '<p style="text-align:center; color: #7f8c8d;">تصویری برای این موجود یافت نشد.</p>';
            }
        }
    }
}

// --- Function to Display Results for Main Search ---
function displayMainSearchResults(data) {
    displayImageForMainSearch(data); 

    resultsContainer.innerHTML = ''; // Clear previous results first

    // Display Wikipedia Summary if available
    if (data.wikipediaSummary) {
        const summaryDiv = document.createElement('div');
        summaryDiv.className = 'wiki-summary';
        // Sanitize summary if needed, but for now, assuming it's safe or will be handled by browser
        summaryDiv.innerHTML = `<p><strong>خلاصه از ویکی‌پدیا:</strong></p><p>${data.wikipediaSummary}</p>`;
        resultsContainer.appendChild(summaryDiv);
    }

    let classificationHtmlOutput = '';
    const hasClassificationDetails = Object.keys(data).some(key => 
        key !== 'imageUrl' && key !== 'searchedName' && key !== 'message' && key !== 'wikipediaSummary' && data[key] !== null && data[key] !== undefined
    );

    if (hasClassificationDetails) {
        classificationHtmlOutput += '<h2>نتایج طبقه‌بندی:</h2>';
        classificationHtmlOutput += '<ul>';

        const displayOrder = [
            { key: 'searchedName', label: 'نام جستجو شده' },
            { key: 'scientificName', label: 'نام علمی' },
            { key: 'kingdom', label: 'سلسله (فرمانرو)' },
            { key: 'phylum', label: 'شاخه' },
            { key: 'class', label: 'رده' },
            { key: 'order', label: 'راسته' },
            { key: 'family', label: 'خانواده' },
            { key: 'genus', label: 'سرده (جنس)' },
            { key: 'species', label: 'گونه' },
            { key: 'rank', label: 'رتبه طبقه‌بندی' },
            { key: 'status', label: 'وضعیت نام' },
            { key: 'matchType', label: 'نوع تطابق GBIF' },
            { key: 'confidence', label: 'درجه اطمینان GBIF (%)' },
        ];

        displayOrder.forEach(item => {
            if (data[item.key] !== undefined && data[item.key] !== null) {
                classificationHtmlOutput += `<li><strong>${item.label}:</strong> ${data[item.key]}</li>`;
            }
        });
        classificationHtmlOutput += '</ul>';
    } else if (data.message && !data.imageUrl && !data.wikipediaSummary) { 
        // Show message only if there's no other content (image, summary, or classification)
        classificationHtmlOutput = `<p class="info-message">${data.message}</p>`;
    } else if (!data.imageUrl && !hasClassificationDetails && !data.wikipediaSummary && !data.message && data.searchedName) { 
         // This case handles when only searchedName is returned, possibly with low confidence match type
        classificationHtmlOutput = `<p class="info-message">اطلاعات دقیقی برای '${data.searchedName}' یافت نشد. لطفاً نام دیگری را امتحان کنید.</p>`;
    }
    
    // Append classification HTML to resultsContainer if it has content
    if (classificationHtmlOutput.trim() !== '') {
        const classificationContentDiv = document.createElement('div');
        classificationContentDiv.innerHTML = classificationHtmlOutput;
        resultsContainer.appendChild(classificationContentDiv);
    }

    // Determine if resultsContainer should be displayed
    if (resultsContainer.innerHTML.trim() !== '') {
        resultsContainer.style.display = 'block';
    } else if (!speciesImageContainer.innerHTML.includes('<img')) { 
        // If no image and no results, hide container
        resultsContainer.style.display = 'none';
    }
    
    errorContainer.style.display = 'none';
}

// --- Function to Show Errors for Main Search ---
function showMainSearchError(message) {
    errorContainer.innerHTML = `<p>${message}</p>`;
    errorContainer.style.display = 'block';
    resultsContainer.innerHTML = ''; 
    resultsContainer.style.display = 'none'; 
    if (speciesImageContainer) speciesImageContainer.innerHTML = '';
    loadingIndicator.style.display = 'none';
}

// --- Function to Show Info Messages for Main Search (in results area) ---
function showMainSearchInfo(message) {
    resultsContainer.innerHTML = `<p class="info-message">${message}</p>`; 
    resultsContainer.style.display = 'block'; 
    errorContainer.innerHTML = '';
    errorContainer.style.display = 'none';
    loadingIndicator.style.display = 'none';
}
