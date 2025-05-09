// URLs for API endpoints
const API_BASE_URL_MAIN_SEARCH = "https://taxonomy-app-ebon.vercel.app/api/handler"; // Your existing main search endpoint
const API_BASE_URL_SUGGESTION = "https://taxonomy-app-ebon.vercel.app/api/suggest_scientific_name"; // New suggestion endpoint

// Elements for Main Search
const speciesNameInput = document.getElementById('speciesNameInput');
const searchButton = document.getElementById('searchButton');
const resultsContainer = document.getElementById('resultsContainer');
const loadingIndicator = document.getElementById('loadingIndicator');
const errorContainer = document.getElementById('errorContainer');
const speciesImageContainer = document.getElementById('speciesImageContainer');

// Elements for Suggestion Feature
const commonNameInput = document.getElementById('commonNameInput');
const suggestNameButton = document.getElementById('suggestNameButton');
const suggestionResultContainer = document.getElementById('suggestionResultContainer');
const suggestionLoadingIndicator = document.getElementById('suggestionLoadingIndicator');

// --- Event Listeners for Main Search ---
if (searchButton) {
    searchButton.addEventListener('click', performMainSearch);
}
if (speciesNameInput) {
    speciesNameInput.addEventListener('keypress', function(event) {
        if (event.key === 'Enter') {
            performMainSearch();
        }
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

// --- Function to Perform Suggestion Search (NCBI) ---
async function performSuggestionSearch() {
    const commonName = commonNameInput.value.trim();
    
    // Clear previous suggestion results and hide indicator
    suggestionResultContainer.innerHTML = '';
    suggestionResultContainer.style.display = 'none';
    suggestionLoadingIndicator.style.display = 'none';

    if (!commonName) {
        showSuggestionMessage("لطفاً نام رایج یک موجود را برای پیشنهاد وارد کنید.", "error");
        return;
    }

    suggestionLoadingIndicator.style.display = 'block'; // Show loading for suggestion

    try {
        const suggestApiUrl = `${API_BASE_URL_SUGGESTION}?common_name=${encodeURIComponent(commonName)}`;
        const response = await fetch(suggestApiUrl);
        const data = await response.json();

        suggestionLoadingIndicator.style.display = 'none'; // Hide loading

        if (!response.ok) {
            // Handle errors from the suggestion API (4xx, 5xx)
            showSuggestionMessage(data.error || data.message || `خطا از سرور پیشنهاد نام (کد: ${response.status})`, "error");
        } else {
            // Success from suggestion API
            if (data.suggested_scientific_name) {
                showSuggestionMessage(
                    `نام علمی پیشنهادی برای '${data.common_name_searched}': <br><strong>${data.suggested_scientific_name}</strong>
                     <button class="copy-to-main-search-btn" data-name="${data.suggested_scientific_name}">کپی به جستجوی اصلی</button>`,
                    "success"
                );
                // Add event listener for the newly created copy button
                const copyBtn = suggestionResultContainer.querySelector('.copy-to-main-search-btn');
                if (copyBtn) {
                    copyBtn.addEventListener('click', function() {
                        speciesNameInput.value = this.dataset.name;
                        speciesNameInput.focus();
                        // Optionally, provide feedback that text was copied
                        showSuggestionMessage(`'${this.dataset.name}' در کادر جستجوی اصلی کپی شد. برای جستجو، دکمه "جستجو" را بزنید.`, "info");
                        // Scroll to main search input for better UX
                        speciesNameInput.scrollIntoView({ behavior: 'smooth', block: 'center' });
                    });
                }
            } else {
                // API returned 200 OK, but no suggestion found (e.g., 404 from NCBI handled by our backend)
                showSuggestionMessage(data.message || `نام علمی برای '${commonName}' پیشنهاد نشد.`, "info");
            }
        }

    } catch (error) {
        suggestionLoadingIndicator.style.display = 'none'; // Hide loading on network error
        showSuggestionMessage(`خطا در برقراری ارتباط با سرور پیشنهاد نام: ${error.message}`, "error");
        console.error("Suggestion Fetch Error:", error);
    }
}

// --- Function to Display Messages for Suggestion Box ---
function showSuggestionMessage(message, type = "info") { // type: "success", "error", "info"
    let cssClass = 'suggestion-message-info'; // Default
    if (type === 'success') {
        cssClass = 'suggested-name-success';
    } else if (type === 'error') {
        cssClass = 'suggestion-message-error';
    }
    // For "info", it will use 'suggestion-message-info'

    suggestionResultContainer.innerHTML = `<p class="${cssClass}">${message}</p>`;
    suggestionResultContainer.style.display = 'block';
}


// --- Function to Perform Main Search (GBIF & Wikipedia) ---
async function performMainSearch() {
    const speciesName = speciesNameInput.value.trim();

    // Clear previous main search results and indicators
    resultsContainer.innerHTML = ''; 
    resultsContainer.style.display = 'none';
    errorContainer.innerHTML = '';   
    errorContainer.style.display = 'none'; 
    if (speciesImageContainer) {
        speciesImageContainer.innerHTML = ''; 
    }
    loadingIndicator.style.display = 'block'; // Show main loading indicator

    if (!speciesName) {
        showMainSearchError("لطفاً نام یک موجود را برای جستجوی اصلی وارد کنید.");
        return;
    }

    try {
        const apiUrl = `${API_BASE_URL_MAIN_SEARCH}?name=${encodeURIComponent(speciesName)}`;
        const response = await fetch(apiUrl);
        const data = await response.json();

        loadingIndicator.style.display = 'none'; // Hide main loading indicator

        if (!response.ok) {
            let errorMessage = data.error || data.message || `خطای ناشناخته از سرور اصلی (کد: ${response.status})`;
             // Special handling for 404s that might still have an image or are just info messages
            if (response.status === 404 && data.message && !data.error) {
                 // If there's an image URL even with a 404 message (e.g., GBIF not found, but Wiki image yes)
                if (data.imageUrl && speciesImageContainer) {
                    displayImageForMainSearch(data); // Display image
                    showMainSearchInfo(data.message); // Show info message in results area
                } else {
                    showMainSearchInfo(data.message); // Show info message in results area, no image
                }
            } else {
                showMainSearchError(errorMessage); 
            }
        } else {
            // Successful response from main search
            if (data.imageUrl && speciesImageContainer) {
                speciesImageContainer.innerHTML = '<div class="loader"></div>'; // Show spinner for image
            } else if (speciesImageContainer) {
                 speciesImageContainer.innerHTML = ''; // Clear if no image URL
            }
            displayMainSearchResults(data);
        }

    } catch (error) {
        showMainSearchError(`خطا در برقراری ارتباط با سرور اصلی: ${error.message}`);
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
                speciesImageContainer.innerHTML = '<p style="text-align:center; color: #721c24;">خطا در بارگذاری تصویر.</p>';
            };
            img.src = data.imageUrl;
            img.alt = data.scientificName || data.searchedName || 'تصویر موجود';
            img.style.maxWidth = "100%";
            img.style.maxHeight = "350px";
            img.style.borderRadius = "8px";
            img.style.display = "block"; // Ensure image is block for auto margins
            img.style.margin = "0 auto 15px auto"; // Center image
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
    displayImageForMainSearch(data); // Handle image display separately

    let htmlOutput = '';
    // Only show title if there are other details apart from image and searchedName
    const hasDetails = Object.keys(data).some(key => 
        key !== 'imageUrl' && key !== 'searchedName' && key !== 'message' && data[key] !== null && data[key] !== undefined
    );

    if (hasDetails) {
        htmlOutput += '<h2>نتایج طبقه‌بندی:</h2>';
        htmlOutput += '<ul>';

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
                htmlOutput += `<li><strong>${item.label}:</strong> ${data[item.key]}</li>`;
            }
        });
        htmlOutput += '</ul>';
    } else if (data.message && !data.imageUrl) { 
        // If only a message (like "GBIF not found") and no image, display it here
        htmlOutput = `<p class="info-message">${data.message}</p>`;
    } else if (!data.imageUrl && !hasDetails && !data.message) {
        // Fallback if no data at all, though this case should be rare if API responds correctly
        htmlOutput = `<p class="info-message">اطلاعاتی برای نمایش یافت نشد.</p>`;
    }


    resultsContainer.innerHTML = htmlOutput;
    if(htmlOutput.trim() !== '') {
        resultsContainer.style.display = 'block';
    } else {
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
    resultsContainer.innerHTML = `<p class="info-message">${message}</p>`; // Uses .info-message style
    resultsContainer.style.display = 'block'; 
    errorContainer.innerHTML = '';
    errorContainer.style.display = 'none';
    // speciesImageContainer might already be handled or cleared by caller
    loadingIndicator.style.display = 'none';
}
