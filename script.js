// URL های API بک‌اند
const API_BASE_URL = "https://taxonomy-app-ebon.vercel.app/api/handler"; // برای جستجوی اصلی
const API_SUGGEST_URL = "https://taxonomy-app-ebon.vercel.app/api/suggest_name"; // برای راهنمای نام

// عناصر HTML اصلی
const speciesNameInput = document.getElementById('speciesNameInput');
const searchButton = document.getElementById('searchButton');
const resultsContainer = document.getElementById('resultsContainer');
const loadingIndicator = document.getElementById('loadingIndicator');
const errorContainer = document.getElementById('errorContainer');
const speciesImageContainer = document.getElementById('speciesImageContainer');

// عناصر HTML راهنما
const commonNameInput = document.getElementById('commonNameInput');
const suggestButton = document.getElementById('suggestButton');
const suggestionResult = document.getElementById('suggestionResult');

// Event Listeners اصلی
searchButton.addEventListener('click', performSearch);
speciesNameInput.addEventListener('keypress', function(event) {
    if (event.key === 'Enter') {
        performSearch();
    }
});

// Event Listeners راهنما
suggestButton.addEventListener('click', getScientificNameSuggestion);
commonNameInput.addEventListener('keypress', function(event) {
    if (event.key === 'Enter') {
        getScientificNameSuggestion();
    }
});

// =============================================
// تابع برای بخش راهنمای نام علمی
// =============================================
async function getScientificNameSuggestion() {
    const commonName = commonNameInput.value.trim();
    if (!suggestionResult) return;

    suggestionResult.innerHTML = '<div class="loader" style="width: 25px; height: 25px; border-width: 3px; border-top-color: #5dade2;"></div>'; // اسپینر کوچکتر با رنگ آبی راهنما
    suggestionResult.classList.remove('error-message'); 

    if (!commonName) {
        suggestionResult.innerHTML = '<span class="error-message">لطفاً نام رایج (انگلیسی) را وارد کنید.</span>';
        return;
    }

    // استفاده از URL کامل بک‌اند
    const suggestApiUrl = `${API_SUGGEST_URL}?query=${encodeURIComponent(commonName)}&lang=en`; 

    try {
        const response = await fetch(suggestApiUrl);
        const data = await response.json();

        if (!response.ok) {
            suggestionResult.innerHTML = `<span class="error-message">${data.message || data.error || 'خطایی در دریافت پیشنهاد رخ داد.'}</span>`;
        } else if (data.scientific_name_suggestion) {
            const suggestion = data.scientific_name_suggestion;
            suggestionResult.innerHTML = `
                <span>پیشنهاد: </span>
                <strong id="suggestedName" style="margin: 0 5px;">${suggestion}</strong>
                <button class="copy-button" onclick="copySuggestionToSearch()" title="کپی در کادر جستجوی اصلی">⬇️ کپی به جستجو</button>
                <button class="copy-button" onclick="copyToClipboard('${suggestion}')" title="کپی در کلیپ‌بورد">📋 کپی</button>
            `;
        } else {
            // این حالت اگر بک‌اند 404 بدهد اتفاق نمی‌افتد، بلکه در !response.ok مدیریت می‌شود
             suggestionResult.innerHTML = `<span class="error-message">پیشنهادی یافت نشد.</span>`;
        }

    } catch (error) {
        console.error("Suggestion Fetch Error:", error);
        suggestionResult.innerHTML = `<span class="error-message">خطا در ارتباط با سرور راهنما.</span>`;
    }
}

// تابع برای کپی کردن پیشنهاد به فیلد جستجوی اصلی
function copySuggestionToSearch() {
    const suggestedNameElement = document.getElementById('suggestedName');
    if (suggestedNameElement && speciesNameInput) {
        speciesNameInput.value = suggestedNameElement.textContent || suggestedNameElement.innerText; 
        speciesNameInput.focus(); 
        // اسکرول به بخش جستجوی اصلی
        const searchBoxElement = document.querySelector('.search-box');
        if (searchBoxElement) {
             searchBoxElement.scrollIntoView({ behavior: 'smooth', block: 'center' });
        }
    }
}

// تابع برای کپی کردن متن به کلیپ‌بورد
function copyToClipboard(text) {
    if (!navigator.clipboard) {
        // Fallback برای مرورگرهای قدیمی‌تر یا محیط‌های ناامن (http)
        try {
            const textArea = document.createElement("textarea");
            textArea.value = text;
            textArea.style.position = "fixed"; // Prevent scrolling to bottom of page in MS Edge.
            textArea.style.opacity = "0";
            document.body.appendChild(textArea);
            textArea.focus();
            textArea.select();
            document.execCommand('copy');
            document.body.removeChild(textArea);
            alert(`"${text}"\nبه کلیپ‌بورد کپی شد! (fallback)`);
        } catch (err) {
             alert('مرورگر شما از کپی خودکار پشتیبانی نمی‌کند. لطفاً دستی کپی کنید.');
        }
        return;
    }
    navigator.clipboard.writeText(text).then(function() {
        alert(`"${text}"\nبه کلیپ‌بورد کپی شد!`);
    }, function(err) {
        console.error('Clipboard copy failed: ', err);
        alert('خطا در کپی متن. ممکن است نیاز به اجازه (permission) باشد.');
    });
}


// =============================================
// توابع برای بخش جستجوی اصلی (مثل قبل)
// =============================================
async function performSearch() {
    const speciesName = speciesNameInput.value.trim();

    resultsContainer.innerHTML = ''; 
    errorContainer.innerHTML = '';   
    errorContainer.style.display = 'none'; 
    if (speciesImageContainer) speciesImageContainer.innerHTML = ''; 
    loadingIndicator.style.display = 'block'; 

    if (!speciesName) {
        showError("لطفاً نام علمی دقیق را وارد کنید.");
        return;
    }

    try {
        // استفاده از URL کامل بک‌اند
        const apiUrl = `${API_BASE_URL}?name=${encodeURIComponent(speciesName)}`;
        const response = await fetch(apiUrl);
        const data = await response.json();

        loadingIndicator.style.display = 'none'; 

        if (!response.ok) {
            let errorMessage = data.error || data.message || `خطای ناشناخته از سرور (کد: ${response.status})`;
            // اینجا دیگر نیازی به تفکیک خطای 404 برای showInfo نیست چون جستجوی اصلی انتظار نام علمی دارد
            showError(errorMessage); 
        } else {
            // قبل از نمایش نتایج، اگر قرار است تصویری لود شود، spinner تصویر را نمایش بده
            if (data.imageUrl && speciesImageContainer) {
                speciesImageContainer.innerHTML = '<div class="loader"></div>'; 
            } else if (speciesImageContainer) {
                 speciesImageContainer.innerHTML = '';
            }
            displayResults(data); 
        }

    } catch (error) {
        showError(`خطا در برقراری ارتباط با سرور: ${error.message}`); 
        console.error("Fetch Error:", error);
    }
}

function displayResults(data) {
    if (speciesImageContainer) {
        if (data.imageUrl) {
            const img = new Image();
            img.onload = function() {
                speciesImageContainer.innerHTML = ''; 
                speciesImageContainer.appendChild(img);
            };
            img.onerror = function() {
                speciesImageContainer.innerHTML = '<p>خطا در بارگذاری تصویر.</p>';
            };
            img.src = data.imageUrl;
            img.alt = data.scientificName || data.searchedName || 'تصویر موجود';
            img.style.maxWidth = "100%";
            img.style.maxHeight = "350px";
            img.style.borderRadius = "8px";
            img.style.boxShadow = "0 4px 8px rgba(0,0,0,0.1)";
        } else {
             if (speciesImageContainer.innerHTML.includes('loader')) {
                 speciesImageContainer.innerHTML = '<p>تصویری برای این موجود یافت نشد.</p>';
            } else if (!data.imageUrl && data.scientificName) { // فقط اگر نتیجه‌ای بود ولی تصویر نداشت
                 speciesImageContainer.innerHTML = '<p>تصویری برای این موجود یافت نشد.</p>';
            }
        }
    }

    let htmlOutput = '<h2>نتایج طبقه‌بندی:</h2>';
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
    
    // اگر پیام خطای GBIF در کنار تصویر آمده بود (حالت خاص در بک‌اند)
    if(data.gbifLookupMessage){
        htmlOutput += `<li><strong>نکته GBIF:</strong> ${data.gbifLookupMessage}</li>`;
    }

    htmlOutput += '</ul>';
    resultsContainer.innerHTML = htmlOutput;
    resultsContainer.style.display = 'block'; 
    errorContainer.style.display = 'none'; 
}

function showError(message) { // برای خطاهای جستجوی اصلی
    errorContainer.innerHTML = `<p>${message}</p>`;
    errorContainer.style.display = 'block';
    resultsContainer.innerHTML = ''; 
    resultsContainer.style.display = 'none'; 
    if (speciesImageContainer) speciesImageContainer.innerHTML = ''; 
    loadingIndicator.style.display = 'none';
}

// تابع showInfo دیگر برای بخش جستجوی اصلی استفاده نمی‌شود
// چون خطای 404 از GBIF را مستقیما به showError می‌فرستیم

/* 
function showInfo(message) { // این تابع ممکن است دیگر لازم نباشد
    resultsContainer.innerHTML = `<p class="info-message">${message}</p>`;
    resultsContainer.style.display = 'block'; 
    errorContainer.innerHTML = '';
    errorContainer.style.display = 'none';
    if (speciesImageContainer) speciesImageContainer.innerHTML = ''; 
    loadingIndicator.style.display = 'none';
} 
*/
