const API_BASE_URL = "https://taxonomy-app-ebon.vercel.app/api/handler"; // آدرس بک‌اند شما

const speciesNameInput = document.getElementById('speciesNameInput');
const searchButton = document.getElementById('searchButton');
const resultsContainer = document.getElementById('resultsContainer');
const loadingIndicator = document.getElementById('loadingIndicator'); // برای لودینگ کلی
const errorContainer = document.getElementById('errorContainer');
const speciesImageContainer = document.getElementById('speciesImageContainer'); // برای تصویر و اسپینر تصویر

searchButton.addEventListener('click', performSearch);
speciesNameInput.addEventListener('keypress', function(event) {
    if (event.key === 'Enter') {
        performSearch();
    }
});

async function performSearch() {
    const speciesName = speciesNameInput.value.trim();

    // پاک‌سازی اولیه
    resultsContainer.innerHTML = ''; 
    errorContainer.innerHTML = '';   
    errorContainer.style.display = 'none'; 
    if (speciesImageContainer) {
        speciesImageContainer.innerHTML = ''; // پاک کردن تصویر یا spinner قبلی
    }
    loadingIndicator.style.display = 'block'; // نمایش نشانگر بارگذاری کلی

    if (!speciesName) {
        showError("لطفاً نام یک موجود را وارد کنید."); // loadingIndicator در showError مخفی می‌شود
        return;
    }

    try {
        const apiUrl = `${API_BASE_URL}?name=${encodeURIComponent(speciesName)}`;
        const response = await fetch(apiUrl);
        const data = await response.json();

        loadingIndicator.style.display = 'none'; // مخفی کردن نشانگر بارگذاری کلی بعد از دریافت پاسخ

        if (!response.ok) {
            let errorMessage = data.error || data.message || `خطای ناشناخته از سرور (کد: ${response.status})`;
            if (response.status === 404 && data.message && !data.error) {
                showInfo(data.message); 
            } else {
                showError(errorMessage); 
            }
        } else {
            // اگر پاسخ موفقیت آمیز بود و قرار است تصویری نمایش داده شود، ابتدا اسپینر تصویر را نمایش بده
            if (data.imageUrl && speciesImageContainer) {
                speciesImageContainer.innerHTML = '<div class="loader"></div>'; // نمایش spinner برای تصویر
            } else if (speciesImageContainer) { // اگر تصویری در کار نیست، بخش تصویر را خالی کن
                 speciesImageContainer.innerHTML = '';
            }
            displayResults(data); // نتایج متنی و تصویر اصلی (اگر بود) در این تابع مدیریت می‌شوند
        }

    } catch (error) {
        showError(`خطا در برقراری ارتباط با سرور: ${error.message}`); // loadingIndicator در showError مخفی می‌شود
        console.error("Fetch Error:", error);
    }
}

function displayResults(data) {
    // بخش تصویر توسط قسمت بالایی در performSearch (نمایش اسپینر) و اینجا (جایگزینی با تصویر) مدیریت می‌شود
    if (speciesImageContainer) {
        if (data.imageUrl) {
            const img = new Image();
            img.onload = function() {
                // وقتی تصویر کاملا لود شد، اسپینر را پاک کرده و تصویر را جایگزین کن
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
            // اگر از بک‌اند imageUrl نیامد، پیام مناسب نمایش بده
            // (این حالت توسط بخش بالایی در performSearch هم پوشش داده میشه اگر imageUrl نباشه)
            if (speciesImageContainer.innerHTML.includes('loader')) { // اگر اسپینر در حال نمایش بود
                 speciesImageContainer.innerHTML = '<p>تصویری برای این موجود یافت نشد.</p>';
            } else if (!data.imageUrl) { // اگر اصلا imageUrl نبود و اسپینری هم نمایش داده نشده بود
                 speciesImageContainer.innerHTML = '<p>تصویری برای این موجود یافت نشد.</p>';
            }
        }
    }

    // نمایش نتایج متنی طبقه‌بندی
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
    if (speciesImageContainer) speciesImageContainer.innerHTML = ''; // پاک کردن تصویر/اسپینر
    loadingIndicator.style.display = 'none';
}

function showInfo(message) {
    resultsContainer.innerHTML = `<p class="info-message">${message}</p>`;
    resultsContainer.style.display = 'block'; 
    errorContainer.innerHTML = '';
    errorContainer.style.display = 'none';
    if (speciesImageContainer) speciesImageContainer.innerHTML = ''; // پاک کردن تصویر/اسپینر
    loadingIndicator.style.display = 'none';
}
