const API_BASE_URL = "https://taxonomy-app-ebon.vercel.app/api/handler"; // آدرس بک‌اند شما

const speciesNameInput = document.getElementById('speciesNameInput');
const searchButton = document.getElementById('searchButton');
const resultsContainer = document.getElementById('resultsContainer');
const loadingIndicator = document.getElementById('loadingIndicator');
const errorContainer = document.getElementById('errorContainer');

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
    loadingIndicator.style.display = 'block';

    if (!speciesName) {
        showError("لطفاً نام یک موجود را وارد کنید.");
        loadingIndicator.style.display = 'none';
        return;
    }

    try {
        const apiUrl = `${API_BASE_URL}?name=${encodeURIComponent(speciesName)}`;
        const response = await fetch(apiUrl);
        loadingIndicator.style.display = 'none';
        const data = await response.json();

        if (!response.ok) {
            // اگر پاسخ سرور خطا بود (مثلا 404, 500)
            // یا اگر خود API ما پیام خطا در بدنه JSON برگردونده (مثل گونه پیدا نشد)
            let errorMessage = data.error || data.message || `خطای ناشناخته از سرور (کد: ${response.status})`;
            showError(errorMessage);
        } else {
            // اگر پاسخ موفقیت آمیز بود و شامل داده‌های طبقه‌بندی بود
            displayResults(data);
        }

    } catch (error) {
        loadingIndicator.style.display = 'none';
        showError(`خطا در برقراری ارتباط با سرور: ${error.message}`);
        console.error("Fetch Error:", error);
    }
}

function displayResults(data) {
    // اگر API ما پیام 'message' برگردانده (یعنی گونه پیدا نشده یا اطمینان کافی نبوده)
    // این حالت توسط response.ok در بالا هم گرفته میشه، اما برای اطمینان بیشتر
    if (data.message && !data.scientificName) { // اگر فقط پیام بود و نه اطلاعات گونه
        showInfo(data.message); // یک تابع جدید برای نمایش پیام‌های اطلاعاتی
        return;
    }

    let htmlOutput = '<h2>نتایج طبقه‌بندی:</h2>';
    // اینجا میتونیم در آینده تصویر رو هم اضافه کنیم
    // htmlOutput += `<div id="speciesImageContainer"></div>`;

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
        // { key: 'usageKey', label: 'GBIF Usage Key' } // این رو معمولا به کاربر نشون نمیدیم
    ];

    displayOrder.forEach(item => {
        if (data[item.key] !== undefined && data[item.key] !== null) {
            htmlOutput += `<li><strong>${item.label}:</strong> ${data[item.key]}</li>`;
        }
    });

    htmlOutput += '</ul>';
    resultsContainer.innerHTML = htmlOutput;

    // اگر usageKey وجود داشت، برای مرحله بعدی (نمایش تصویر) آماده میشیم
    // if (data.usageKey) {
    //     fetchAndDisplayImage(data.usageKey);
    // }
}

function showError(message) {
    errorContainer.innerHTML = `<p>${message}</p>`;
    errorContainer.style.display = 'block';
    resultsContainer.innerHTML = '';
    loadingIndicator.style.display = 'none'; // مطمئن بشیم که لودینگ هم مخفیه
}

function showInfo(message) { // تابع جدید برای پیام‌های اطلاعاتی
    resultsContainer.innerHTML = `<p class="info-message">${message}</p>`; // یک کلاس برای استایل متفاوت
    errorContainer.innerHTML = '';
    errorContainer.style.display = 'none';
    loadingIndicator.style.display = 'none';
}

// تابع برای گرفتن و نمایش تصویر (در مرحله بعدی پیاده‌سازی میشه)
// async function fetchAndDisplayImage(usageKey) {
//     const imageContainer = document.getElementById('speciesImageContainer');
//     if (!imageContainer) return;
//     imageContainer.innerHTML = '<p>در حال بارگذاری تصویر...</p>';
//     // ... کد مربوط به API تصاویر GBIF ...
// }
