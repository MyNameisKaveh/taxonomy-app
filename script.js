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
const autocompleteSuggestionsContainer = document.getElementById('autocompleteSuggestions');

// Elements for Suggestion Feature
const commonNameInput = document.getElementById('commonNameInput');
const suggestNameButton = document.getElementById('suggestNameButton');
const suggestionResultContainer = document.getElementById('suggestionResultContainer');
const suggestionLoadingIndicator = document.getElementById('suggestionLoadingIndicator');

// Elements for Search History
const searchHistoryList = document.getElementById('searchHistoryList');
const clearHistoryButton = document.getElementById('clearHistoryButton');

// Constants for Search History
const MAX_HISTORY_ITEMS = 10;
const HISTORY_STORAGE_KEY = 'taxonomySearchHistory';

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
            if (autocompleteSuggestionsContainer.style.display === 'block' && autocompleteSuggestionsContainer.firstChild) {
                autocompleteSuggestionsContainer.innerHTML = '';
                autocompleteSuggestionsContainer.style.display = 'none';
            }
            performMainSearch();
        }
    });
    speciesNameInput.addEventListener('input', debounce(handleAutocompleteInput, 300));
    speciesNameInput.addEventListener('blur', function() {
        setTimeout(() => {
            if (!autocompleteSuggestionsContainer.matches(':hover')) {
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

// --- Search History Logic ---
function loadSearchHistory() {
    const historyJson = localStorage.getItem(HISTORY_STORAGE_KEY);
    const historyArray = historyJson ? JSON.parse(historyJson) : [];
    renderSearchHistory(historyArray);
}

function renderSearchHistory(historyArray) {
    searchHistoryList.innerHTML = ''; // Clear existing items

    if (!historyArray || historyArray.length === 0) {
        const emptyMsgLi = document.createElement('li');
        emptyMsgLi.textContent = 'تاریخچه جستجو خالی است.';
        emptyMsgLi.style.fontStyle = 'italic';
        emptyMsgLi.style.color = '#6c757d'; // Using a muted color
        emptyMsgLi.style.cursor = 'default';
        searchHistoryList.appendChild(emptyMsgLi);
        clearHistoryButton.style.display = 'none'; // Hide button if history is empty
        return;
    }

    // Slice to ensure we only render MAX_HISTORY_ITEMS, taking the newest ones (from the start of stored array)
    const itemsToRender = historyArray.slice(0, MAX_HISTORY_ITEMS);

    itemsToRender.forEach(term => {
        const li = document.createElement('li');
        li.textContent = term;
        li.dataset.searchTerm = term;
        li.addEventListener('click', function() {
            speciesNameInput.value = this.dataset.searchTerm;
            // Optionally, scroll to top or to the input field
            speciesNameInput.scrollIntoView({ behavior: 'smooth', block: 'center' });
            performMainSearch();
            // Moving term to top is handled by addSearchToHistory called in performMainSearch
        });
        searchHistoryList.appendChild(li); // Appending, so newest (from unshift) will be at the top visually if CSS is default
                                          // If we prepend here, the order from localStorage (newest first) is maintained.
                                          // Let's stick to addSearchToHistory unshifting, and here we iterate and append.
                                          // This means the list will show newest at top if array is newest first.
    });
    
    // Correcting display order if using appendChild: iterate reversed or prepend.
    // Since addSearchToHistory uses unshift (adds to beginning), the historyArray is newest first.
    // So, when appending, the oldest of the MAX_HISTORY_ITEMS will be at top.
    // To show newest at top with appendChild, we'd iterate reversed.
    // OR, using prepend here:
    searchHistoryList.innerHTML = ''; // Clear again before prepending
    itemsToRender.forEach(term => {
        const li = document.createElement('li');
        li.textContent = term;
        li.dataset.searchTerm = term;
        li.addEventListener('click', function() {
            speciesNameInput.value = this.dataset.searchTerm;
            speciesNameInput.scrollIntoView({ behavior: 'smooth', block: 'center' });
            performMainSearch();
        });
        searchHistoryList.prepend(li); // Prepend to show newest items at the top
    });


    clearHistoryButton.style.display = 'block'; // Show button if history is not empty
}

function addSearchToHistory(searchTerm) {
    if (!searchTerm) return; // Do not add empty searches
    let historyArray = JSON.parse(localStorage.getItem(HISTORY_STORAGE_KEY) || '[]');
    
    // Remove term if it already exists to move it to the top (most recent)
    historyArray = historyArray.filter(item => item.toLowerCase() !== searchTerm.toLowerCase());
    
    historyArray.unshift(searchTerm); // Add new term to the beginning (most recent)
    historyArray = historyArray.slice(0, MAX_HISTORY_ITEMS); // Limit history size
    
    localStorage.setItem(HISTORY_STORAGE_KEY, JSON.stringify(historyArray));
    renderSearchHistory(historyArray);
}

if (clearHistoryButton) {
    clearHistoryButton.addEventListener('click', function() {
        localStorage.removeItem(HISTORY_STORAGE_KEY);
        renderSearchHistory([]); // Render an empty list
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
    autocompleteSuggestionsContainer.innerHTML = '<div>درحال جستجو...</div>';
    autocompleteSuggestionsContainer.style.display = 'block';
    try {
        const response = await fetch(`${GBIF_SUGGEST_API_URL}?q=${encodeURIComponent(trimmedQuery)}&limit=7`);
        if (!response.ok) throw new Error(`GBIF API error: ${response.status}`);
        const suggestions = await response.json();
        autocompleteSuggestionsContainer.innerHTML = '';
        if (suggestions && suggestions.length > 0) {
            suggestions.forEach(suggestion => {
                const itemDiv = document.createElement('div');
                const nameToDisplay = suggestion.scientificName || suggestion.canonicalName || suggestion.vernacularName;
                let displayText = nameToDisplay;
                if (suggestion.scientificName && suggestion.vernacularName) displayText = `<em>${suggestion.scientificName}</em> (${suggestion.vernacularName})`;
                else if (suggestion.canonicalName && suggestion.vernacularName) displayText = `<em>${suggestion.canonicalName}</em> (${suggestion.vernacularName})`;
                else if (suggestion.scientificName) displayText = `<em>${suggestion.scientificName}</em>`;
                else if (suggestion.canonicalName) displayText = `<em>${suggestion.canonicalName}</em>`;

                if (nameToDisplay) {
                    itemDiv.innerHTML = displayText;
                    itemDiv.addEventListener('click', function() {
                        speciesNameInput.value = suggestion.scientificName || suggestion.canonicalName;
                        autocompleteSuggestionsContainer.innerHTML = '';
                        autocompleteSuggestionsContainer.style.display = 'none';
                        performMainSearch();
                    });
                    autocompleteSuggestionsContainer.appendChild(itemDiv);
                }
            });
            if (autocompleteSuggestionsContainer.children.length > 0) autocompleteSuggestionsContainer.style.display = 'block';
            else autocompleteSuggestionsContainer.innerHTML = '<div>نتیجه‌ای یافت نشد.</div>';
        } else {
            autocompleteSuggestionsContainer.innerHTML = '<div>نتیجه‌ای یافت نشد.</div>';
        }
    } catch (error) {
        console.error("Autocomplete fetch error:", error);
        autocompleteSuggestionsContainer.innerHTML = '<div>خطا در دریافت پیشنهادات.</div>';
    }
    if (autocompleteSuggestionsContainer.children.length === 0) autocompleteSuggestionsContainer.style.display = 'none';
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
            showSuggestionMessage(data.error || data.message || `مشکلی در ارتباط با سرویس پیشنهاد نام رخ داد (کد: ${response.status}). لطفاً دوباره تلاش کنید.`, "error");
        } else {
            if (data.suggested_scientific_name) {
                showSuggestionMessage(`نام علمی پیشنهادی برای '${data.common_name_searched}': <br><strong>${data.suggested_scientific_name}</strong> <button class="copy-to-main-search-btn" data-name="${data.suggested_scientific_name}">استفاده از این نام</button>`, "success");
                const copyBtn = suggestionResultContainer.querySelector('.copy-to-main-search-btn');
                if (copyBtn) {
                    copyBtn.addEventListener('click', function() {
                        speciesNameInput.value = this.dataset.name;
                        speciesNameInput.focus();
                        showSuggestionMessage(`نام علمی '${this.dataset.name}' به کادر جستجوی اصلی منتقل شد. اکنون می‌توانید دکمه "جستجو" را بزنید.`, "info");
                        speciesNameInput.scrollIntoView({ behavior: 'smooth', block: 'center' });
                    });
                }
            } else {
                showSuggestionMessage(data.message || `متاسفانه، نام علمی دقیقی برای '${commonName}' از طریق NCBI یافت نشد. لطفاً نام دیگری را امتحان کنید یا از نام علمی شناخته شده استفاده نمایید.`, "info");
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
    if (type === 'success') cssClass = 'suggested-name-success';
    else if (type === 'error') cssClass = 'suggestion-message-error';
    suggestionResultContainer.innerHTML = `<p class="${cssClass}">${message}</p>`;
    suggestionResultContainer.style.display = 'block';
}

// --- Function to Perform Main Search (GBIF & Wikipedia) ---
async function performMainSearch() {
    const speciesName = speciesNameInput.value.trim();

    if (autocompleteSuggestionsContainer) {
        autocompleteSuggestionsContainer.innerHTML = '';
        autocompleteSuggestionsContainer.style.display = 'none';
    }

    if (!speciesName) {
        showMainSearchError("لطفاً نام یک موجود را برای جستجوی اطلاعات وارد کنید.");
        return; 
    }
    
    addSearchToHistory(speciesName); // Add to history

    resultsContainer.innerHTML = '';
    resultsContainer.style.display = 'none';
    errorContainer.innerHTML = '';
    errorContainer.style.display = 'none';
    if (speciesImageContainer) speciesImageContainer.innerHTML = '';
    loadingIndicator.style.display = 'block';

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
            if (data.imageUrl && speciesImageContainer) speciesImageContainer.innerHTML = '<div class="loader"></div>';
            else if (speciesImageContainer) speciesImageContainer.innerHTML = '';
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
            if (speciesImageContainer.innerHTML.includes('loader')) speciesImageContainer.innerHTML = '<p style="text-align:center; color: #7f8c8d;">تصویری برای این موجود یافت نشد.</p>';
            else speciesImageContainer.innerHTML = '<p style="text-align:center; color: #7f8c8d;">تصویری برای این موجود یافت نشد.</p>';
        }
    }
}

// --- Function to Display Results for Main Search ---
function displayMainSearchResults(data) {
    displayImageForMainSearch(data);
    resultsContainer.innerHTML = '';
    if (data.wikipediaSummary) {
        const summaryDiv = document.createElement('div');
        summaryDiv.className = 'wiki-summary';
        summaryDiv.innerHTML = `<p><strong>خلاصه از ویکی‌پدیا:</strong></p><p>${data.wikipediaSummary}</p>`;
        resultsContainer.appendChild(summaryDiv);
    }
    let classificationHtmlOutput = '';
    const hasClassificationDetails = Object.keys(data).some(key => key !== 'imageUrl' && key !== 'searchedName' && key !== 'message' && key !== 'wikipediaSummary' && data[key] !== null && data[key] !== undefined);
    if (hasClassificationDetails) {
        classificationHtmlOutput += '<h2>نتایج طبقه‌بندی:</h2><ul>';
        const displayOrder = [
            { key: 'searchedName', label: 'نام جستجو شده' }, { key: 'scientificName', label: 'نام علمی' },
            { key: 'kingdom', label: 'سلسله (فرمانرو)' }, { key: 'phylum', label: 'شاخه' },
            { key: 'class', label: 'رده' }, { key: 'order', label: 'راسته' },
            { key: 'family', label: 'خانواده' }, { key: 'genus', label: 'سرده (جنس)' },
            { key: 'species', label: 'گونه' }, { key: 'rank', label: 'رتبه طبقه‌بندی' },
            { key: 'status', label: 'وضعیت نام' }, { key: 'matchType', label: 'نوع تطابق GBIF' },
            { key: 'confidence', label: 'درجه اطمینان GBIF (%)' },
        ];
        displayOrder.forEach(item => {
            if (data[item.key] !== undefined && data[item.key] !== null) {
                classificationHtmlOutput += `<li><strong>${item.label}:</strong> ${data[item.key]}</li>`;
            }
        });
        classificationHtmlOutput += '</ul>';
    } else if (data.message && !data.imageUrl && !data.wikipediaSummary) {
        classificationHtmlOutput = `<p class="info-message">${data.message}</p>`;
    } else if (!data.imageUrl && !hasClassificationDetails && !data.wikipediaSummary && !data.message && data.searchedName) {
        classificationHtmlOutput = `<p class="info-message">اطلاعات دقیقی برای '${data.searchedName}' یافت نشد. لطفاً نام دیگری را امتحان کنید.</p>`;
    }
    if (classificationHtmlOutput.trim() !== '') {
        const classificationContentDiv = document.createElement('div');
        classificationContentDiv.innerHTML = classificationHtmlOutput;
        resultsContainer.appendChild(classificationContentDiv);
    }
    if (resultsContainer.innerHTML.trim() !== '') resultsContainer.style.display = 'block';
    else if (!speciesImageContainer.innerHTML.includes('<img')) resultsContainer.style.display = 'none';
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

// --- Initializations ---
loadSearchHistory(); // Load search history on page load
