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

    resultsContainer.innerHTML = ''; // پاک کردن نتایج قبلی
    errorContainer.innerHTML = '';   // پاک کردن خطاهای قبلی
    errorContainer.style.display = 'none'; // مخفی کردن کادر خطا
    loadingIndicator.style.display = 'block'; // نمایش نشانگر بارگذاری

    if (!speciesName) {
        showError("لطفاً نام یک موجود را وارد کنید.");
        // loadingIndicator.style.display = 'none'; // این خط در showError انجام میشه
        return;
    }

    try {
        const apiUrl = `${API_BASE_URL}?name=${encodeURIComponent(speciesName)}`;
        const response = await fetch(apiUrl);
        // loadingIndicator.style.display = 'none'; // این خط هم در showError یا بعد از گرفتن دیتا انجام میشه
        const data = await response.json();

        if (!response.ok) {
            // اگر پاسخ سرور خطا بود (مثلا 404, 500)
            let errorMessage = data.error || data.message || `خطای ناشناخته از سرور (کد: ${response.status})`;
            
            // چک می‌کنیم که آیا این "خطا" در واقع پیام "گونه پیدا نشد" است یا نه
            // بک‌اند ما برای "گونه پیدا نشد" کد 404 و یک فیلد 'message' برمی‌گردونه
            if (response.status === 404 && data.message && !data.error) {
                loadingIndicator.style.display = 'none'; // مخفی کردن لودینگ قبل از نمایش اطلاعات
                showInfo(data.message); // اگر فقط پیام "گونه پیدا نشد" بود، با showInfo نمایش بده
            } else {
                showError(errorMessage); // در غیر این صورت، به عنوان خطا نمایش بده (لودینگ در showError مخفی میشه)
            }
        } else {
            // اگر همه چیز موفقیت‌آمیز بود و شامل داده‌های طبقه‌بندی بود
            loadingIndicator.style.display = 'none'; // مخفی کردن لودینگ قبل از نمایش نتایج
            displayResults(data);
        }

    } catch (error) {
        // loadingIndicator.style.display = 'none'; // این خط در showError انجام میشه
        showError(`خطا در برقراری ارتباط با سرور: ${error.message}`);
        console.error("Fetch Error:", error);
    }
}

function displayResults(data) {
    // این شرط الان در performSearch مدیریت میشه، اما برای اطمینان بیشتر می‌تونه بمونه
    // if (data.message && !data.scientificName) {
    //     showInfo(data.message);
    //     return;
    // }

    let htmlOutput = '<h2>نتایج طبقه‌بندی:</h2>';
    // اینجا در آینده تصویر رو اضافه می‌کنیم
    // htmlOutput += `<div id="speciesImageContainer" style="text-align: center; margin-bottom: 15px;"></div>`;

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
    resultsContainer.style.display = 'block'; // مطمئن بشیم که کانتینر نتایج نمایش داده میشه
    errorContainer.style.display = 'none'; // مخفی کردن کادر خطا در صورت موفقیت

    // اگر usageKey وجود داشت، برای مرحله بعدی (نمایش تصویر) آماده میشیم
    // if (data.usageKey) {
    //     fetchAndDisplayImage(data.usageKey);
    // }
}

function showError(message) {
    errorContainer.innerHTML = `<p>${message}</p>`;
    errorContainer.style.display = 'block';
    resultsContainer.innerHTML = ''; // پاک کردن نتایج قبلی اگر خطایی رخ داد
    resultsContainer.style.display = 'none'; // مخفی کردن کانتینر نتایج
    loadingIndicator.style.display = 'none';
}

function showInfo(message) {
    resultsContainer.innerHTML = `<p class="info-message">${message}</p>`;
    resultsContainer.style.display = 'block'; // نمایش کانتینر نتایج برای پیام اطلاعاتی
    errorContainer.innerHTML = '';
    errorContainer.style.display = 'none';
    loadingIndicator.style.display = 'none';
}

// تابع برای گرفتن و نمایش تصویر (در مرحله بعدی پیاده‌سازی میشه)
// async function fetchAndDisplayImage(usageKey) {
//     const imageContainer = document.getElementById('speciesImageContainer');
//     if (!imageContainer) return; // اگر در HTML این بخش رو اضافه نکرده باشیم
//     imageContainer.innerHTML = '<p>در حال بارگذاری تصویر...</p>';
//     try {
//         // آدرس API تصاویر GBIF: https://api.gbif.org/v1/species/{usageKey}/media
//         const imageApiUrl = `https://api.gbif.org/v1/species/${usageKey}/media?limit=1`; // فقط یک تصویر
//         const response = await fetch(imageApiUrl);
//         if (!response.ok) {
//             imageContainer.innerHTML = '<p>تصویری برای این موجود یافت نشد.</p>';
//             return;
//         }
//         const mediaData = await response.json();
//         if (mediaData.results && mediaData.results.length > 0 && mediaData.results[0].identifier) {
//             const imageUrl = mediaData.results[0].identifier;
//             // برخی URLها کامل نیستند، باید چک کنیم
//             const finalImageUrl = imageUrl.startsWith('http') ? imageUrl : `https:${imageUrl}`; // اگر http نداشت، https رو اضافه کنیم (یا از خود GBIF بپرسیم)
//             imageContainer.innerHTML = `<img src="${finalImageUrl}" alt="تصویر ${speciesNameInput.value}" style="max-width: 100%; max-height: 300px; border-radius: 4px;">`;
//         } else {
//             imageContainer.innerHTML = '<p>تصویری برای این موجود یافت نشد.</p>';
//         }
//     } catch (error) {
//         console.error("Error fetching image:", error);
//         imageContainer.innerHTML = '<p>خطا در بارگذاری تصویر.</p>';
//     }
// }
