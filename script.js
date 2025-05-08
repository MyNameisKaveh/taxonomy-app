const API_BASE_URL = "https://taxonomy-app-ebon.vercel.app/api/handler";

const speciesNameInput = document.getElementById('speciesNameInput');
const searchButton = document.getElementById('searchButton');
const resultsContainer = document.getElementById('resultsContainer');
const loadingIndicator = document.getElementById('loadingIndicator');
const errorContainer = document.getElementById('errorContainer');
const speciesImageContainer = document.getElementById('speciesImageContainer'); // عنصر جدید برای تصویر

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
    speciesImageContainer.innerHTML = ''; // پاک کردن تصویر قبلی
    loadingIndicator.style.display = 'block';

    if (!speciesName) {
        showError("لطفاً نام یک موجود را وارد کنید.");
        return;
    }

    try {
        const apiUrl = `${API_BASE_URL}?name=${encodeURIComponent(speciesName)}`;
        const response = await fetch(apiUrl);
        const data = await response.json();

        loadingIndicator.style.display = 'none'; // مخفی کردن لودینگ بعد از دریافت پاسخ

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

    // اگر usageKey وجود داشت، تصویر را بارگذاری کن
    if (data.usageKey) {
        fetchAndDisplayImage(data.usageKey, data.scientificName || speciesNameInput.value);
    } else {
        speciesImageContainer.innerHTML = ''; // اگر کلیدی نبود، بخش تصویر را خالی کن
    }
}

async function fetchAndDisplayImage(usageKey, altText = "تصویر موجود") {
    if (!speciesImageContainer) return;
    speciesImageContainer.innerHTML = '<p>در حال بارگذاری تصویر...</p>';

    try {
        // آدرس API تصاویر GBIF: https://api.gbif.org/v1/species/{usageKey}/media
        // limit=1 برای گرفتن فقط اولین تصویر
        // mediaType=StillImage برای اینکه فقط تصاویر ثابت رو بگیریم (نه ویدیو یا صدا)
        const imageApiUrl = `https://api.gbif.org/v1/species/${usageKey}/media?limit=1&mediaType=StillImage`;
        const response = await fetch(imageApiUrl);

        if (!response.ok) {
            // اگر API تصاویر خطایی برگرداند یا نتیجه‌ای نداشت
            console.warn(`GBIF Media API for key ${usageKey} returned status ${response.status}`);
            speciesImageContainer.innerHTML = '<p>تصویری برای این موجود در GBIF یافت نشد.</p>';
            return;
        }

        const mediaData = await response.json();

        if (mediaData.results && mediaData.results.length > 0 && mediaData.results[0].identifier) {
            let imageUrl = mediaData.results[0].identifier;
            
            // برخی URL های تصاویر در GBIF ممکن است کامل نباشند یا نیاز به تغییر اندازه داشته باشند
            // اینجا یک URL تصویر کوچک (thumbnail) یا متوسط رو انتخاب می‌کنیم اگر موجود باشه
            // یا اگر URL کامل نبود، https به آن اضافه می‌کنیم
            // این بخش ممکنه نیاز به تنظیم دقیق‌تر بر اساس پاسخ GBIF داشته باشه

            if (mediaData.results[0].references) { // گاهی URL اصلی در references است
                imageUrl = mediaData.results[0].references;
            }

            // اگر URL با // شروع میشه، https: رو بهش اضافه کن
            if (imageUrl.startsWith("//")) {
                imageUrl = "https:" + imageUrl;
            }
            // اگر با http یا https شروع نمیشد، و یک آدرس نسبی بود، ممکنه کار نکنه.
            // GBIF معمولا URL های کامل برمیگردونه.

            speciesImageContainer.innerHTML = `<img src="${imageUrl}" alt="${altText}" style="max-width: 100%; max-height: 350px; border-radius: 8px; box-shadow: 0 4px 8px rgba(0,0,0,0.1);">`;
        } else {
            speciesImageContainer.innerHTML = '<p>تصویری برای این موجود در GBIF یافت نشد.</p>';
        }
    } catch (error) {
        console.error("Error fetching or displaying image:", error);
        speciesImageContainer.innerHTML = '<p>خطا در بارگذاری تصویر.</p>';
    }
}

function showError(message) {
    errorContainer.innerHTML = `<p>${message}</p>`;
    errorContainer.style.display = 'block';
    resultsContainer.innerHTML = '';
    resultsContainer.style.display = 'none';
    speciesImageContainer.innerHTML = ''; // پاک کردن تصویر در صورت خطا
    loadingIndicator.style.display = 'none';
}

function showInfo(message) {
    resultsContainer.innerHTML = `<p class="info-message">${message}</p>`;
    resultsContainer.style.display = 'block';
    errorContainer.innerHTML = '';
    errorContainer.style.display = 'none';
    speciesImageContainer.innerHTML = ''; // پاک کردن تصویر برای پیام اطلاعاتی
    loadingIndicator.style.display = 'none';
}
