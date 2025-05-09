// URLs for API endpoints
const API_BASE_URL_MAIN_SEARCH = "https://taxonomy-app-ebon.vercel.app/api/handler"; 
const API_BASE_URL_SUGGESTION = "https://taxonomy-app-ebon.vercel.app/api/suggest_scientific_name"; 

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
            ); // متن خطا بهبود یافت
        } else {
            if (data.suggested_scientific_name) {
                showSuggestionMessage(
                    `نام علمی پیشنهادی برای '${data.common_name_searched}': <br><strong>${data.suggested_scientific_name}</strong>
                     <button class="copy-to-main-search-btn" data-name="${data.suggested_scientific_name}">استفاده از این نام</button>`, // متن دکمه بهبود یافت
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
                        ); // متن پیام بهبود یافت
                        speciesNameInput.scrollIntoView({ behavior: 'smooth', block: 'center' });
                    });
                }
            } else {
                showSuggestionMessage(
                    data.message || `متاسفانه، نام علمی دقیقی برای '${commonName}' از طریق NCBI یافت نشد. لطفاً نام دیگری را امتحان کنید یا از نام علمی شناخته شده استفاده نمایید.`, 
                    "info"
                ); // متن پیام بهبود یافت
            }
        }

    } catch (error) {
        suggestionLoadingIndicator.style.display = 'none'; 
        showSuggestionMessage(`مشکلی در ارتباط با سرویس پیشنهاد نام رخ داد. لطفاً اتصال اینترنت خود را بررسی کرده و دوباره تلاش کنید. (پیام سیستم: ${error.message})`, "error"); // متن خطا بهبود یافت
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

    resultsContainer.innerHTML = ''; 
    resultsContainer.style.display = 'none';
    errorContainer.innerHTML = '';   
    errorContainer.style.display = 'none'; 
    if (speciesImageContainer) {
        speciesImageContainer.innerHTML = ''; 
    }
    loadingIndicator.style.display = 'block'; 

    if (!speciesName) {
        showMainSearchError("لطفاً نام یک موجود را برای جستجوی اطلاعات وارد کنید."); // متن بهبود یافت
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
                    showMainSearchInfo(data.message || `اطلاعات طبقه‌بندی برای '${speciesName}' در GBIF یافت نشد یا دقت کافی نداشت. اگر نام رایج وارد کرده‌اید، سعی کنید نام علمی آن را پیدا و جستجو کنید.`); // متن بهبود یافت
                } else {
                    showMainSearchInfo(data.message || `اطلاعات طبقه‌بندی برای '${speciesName}' در GBIF یافت نشد یا دقت کافی نداشت. اگر نام رایج وارد کرده‌اید، سعی کنید نام علمی آن را پیدا و جستجو کنید.`); // متن بهبود یافت
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
        showMainSearchError(`مشکلی در ارتباط با سرور اصلی رخ داد. لطفاً اتصال اینترنت خود را بررسی کرده و دوباره تلاش کنید. (پیام سیستم: ${error.message})`); // متن خطا بهبود یافت
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
                speciesImageContainer.innerHTML = '<p style="text-align:center; color: #721c24;">خطا در بارگذاری تصویر. ممکن است آدرس تصویر نامعتبر باشد.</p>'; // متن بهبود یافت
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

    let htmlOutput = '';
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
        htmlOutput = `<p class="info-message">${data.message}</p>`;
    } else if (!data.imageUrl && !hasDetails && !data.message && Object.keys(data).length <=2) { // <=2 برای searchedName و matchType/confidence
        htmlOutput = `<p class="info-message">اطلاعات دقیقی برای '${data.searchedName || 'این موجود'}' یافت نشد. لطفاً نام دیگری را امتحان کنید.</p>`; // متن بهبود یافت
    }


    resultsContainer.innerHTML = htmlOutput;
    if(htmlOutput.trim() !== '') {
        resultsContainer.style.display = 'block';
    } else if (!speciesImageContainer.innerHTML.includes('<img')) { // اگر تصویری هم نبود، مخفی کن
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
