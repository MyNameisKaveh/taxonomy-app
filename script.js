const API_BASE_URL = "https://taxonomy-app-ebon.vercel.app/api/handler"; // آدرس بک‌اند شما

const speciesNameInput = document.getElementById('speciesNameInput');
const searchButton = document.getElementById('searchButton');
const resultsContainer = document.getElementById('resultsContainer');
const loadingIndicator = document.getElementById('loadingIndicator');
const errorContainer = document.getElementById('errorContainer');
const speciesImageContainer = document.getElementById('speciesImageContainer'); // عنصر برای تصویر

searchButton.addEventListener('click', performSearch);
speciesNameInput.addEventListener('keypress', function(event) {
    if (event.key === 'Enter') {
        performSearch();
    }
});

async function performSearch() {
    const speciesName = speciesNameInput.value.trim();

    resultsContainer.innerHTML = ''; 
    errorContainer.innerHTML = '';   
    errorContainer.style.display = 'none'; 
    if (speciesImageContainer) speciesImageContainer.innerHTML = ''; // پاک کردن تصویر قبلی
    loadingIndicator.style.display = 'block'; 

    if (!speciesName) {
        showError("لطفاً نام یک موجود را وارد کنید.");
        return;
    }

    try {
        const apiUrl = `${API_BASE_URL}?name=${encodeURIComponent(speciesName)}`;
        const response = await fetch(apiUrl);
        const data = await response.json();

        loadingIndicator.style.display = 'none'; 

        if (!response.ok) {
            let errorMessage = data.error || data.message || `خطای ناشناخته از سرور (کد: ${response.status})`;
            if (response.status === 404 && data.message && !data.error) {
                showInfo(data.message); 
            } else {
                showError(errorMessage); 
            }
        } else {
            displayResults(data);
        }

    } catch (error) {
        showError(`خطا در برقراری ارتباط با سرور: ${error.message}`);
        console.error("Fetch Error:", error);
    }
}

function displayResults(data) {
    // نمایش تصویر اگر imageUrl در پاسخ بک‌اند وجود داشت
    if (speciesImageContainer) { // ابتدا مطمئن شویم عنصر تصویر در DOM وجود دارد
        if (data.imageUrl) {
            speciesImageContainer.innerHTML = `<img src="${data.imageUrl}" alt="${data.scientificName || data.searchedName || 'تصویر موجود'}" style="max-width: 100%; max-height: 350px; border-radius: 8px; box-shadow: 0 4px 8px rgba(0,0,0,0.1);">`;
        } else {
            speciesImageContainer.innerHTML = '<p>تصویری برای این موجود یافت نشد.</p>';
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

    htmlOutput += '</ul>';
    resultsContainer.innerHTML = htmlOutput;
    resultsContainer.style.display = 'block'; 
    errorContainer.style.display = 'none'; 
}

function showError(message) {
    errorContainer.innerHTML = `<p>${message}</p>`;
    errorContainer.style.display = 'block';
    resultsContainer.innerHTML = ''; 
    resultsContainer.style.display = 'none'; 
    if (speciesImageContainer) speciesImageContainer.innerHTML = ''; // پاک کردن تصویر در صورت خطا
    loadingIndicator.style.display = 'none';
}

function showInfo(message) {
    // اگر میخواهید پیام اطلاعاتی (مثل گونه پیدا نشد) در بخش نتایج نمایش داده شود
    // و با استایل متفاوتی (مثلا آبی)
    resultsContainer.innerHTML = `<p class="info-message">${message}</p>`;
    resultsContainer.style.display = 'block'; 
    errorContainer.innerHTML = '';
    errorContainer.style.display = 'none';
    if (speciesImageContainer) speciesImageContainer.innerHTML = ''; // پاک کردن تصویر برای پیام اطلاعاتی
    loadingIndicator.style.display = 'none';
}
