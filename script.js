// آدرس API بک‌اند شما در Vercel
const API_BASE_URL = "https://taxonomy-app-ebon.vercel.app/api/handler"; // !!! آدرس خودت رو اینجا جایگزین کن !!!

// گرفتن عناصر HTML از صفحه
const speciesNameInput = document.getElementById('speciesNameInput');
const searchButton = document.getElementById('searchButton');
const resultsContainer = document.getElementById('resultsContainer');
const loadingIndicator = document.getElementById('loadingIndicator');
const errorContainer = document.getElementById('errorContainer');

// اضافه کردن event listener به دکمه جستجو
searchButton.addEventListener('click', performSearch);

// اضافه کردن event listener به فیلد ورودی برای جستجو با زدن Enter
speciesNameInput.addEventListener('keypress', function(event) {
    if (event.key === 'Enter') {
        performSearch();
    }
});

// تابع اصلی برای انجام جستجو
async function performSearch() {
    const speciesName = speciesNameInput.value.trim(); // گرفتن مقدار ورودی و حذف فاصله‌های اضافی ابتدا و انتها

    // پاک کردن نتایج و خطاهای قبلی و نمایش نشانگر بارگذاری
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
        // ساخت URL کامل برای درخواست به API
        // استفاده از encodeURIComponent برای مدیریت کاراکترهای خاص مثل فاصله در نام گونه
        const apiUrl = `${API_BASE_URL}?name=${encodeURIComponent(speciesName)}`;

        const response = await fetch(apiUrl, {
            method: 'GET' // متد GET برای API ما کافی است
        });

        loadingIndicator.style.display = 'none'; // مخفی کردن نشانگر بارگذاری

        const data = await response.json(); // تبدیل پاسخ به JSON

        if (!response.ok) {
            // اگر پاسخ سرور خطا بود (مثلا 404, 500)
            showError(`خطا از سرور: ${data.error || response.statusText}`);
        } else {
            // اگر همه چیز موفقیت‌آمیز بود
            displayResults(data);
        }

    } catch (error) {
        // برای خطاهای شبکه یا خطاهای دیگر در زمان اجرای fetch
        loadingIndicator.style.display = 'none';
        showError(`خطا در برقراری ارتباط: ${error.message}`);
        console.error("Fetch Error:", error);
    }
}

// تابع برای نمایش نتایج در صفحه
function displayResults(data) {
    if (data.error) { // اگر API خودمون خطایی برگردونده بود (مثلا گونه پیدا نشد)
        showError(data.error);
        return;
    }

    let htmlOutput = '<h2>نتایج طبقه‌بندی:</h2>';
    htmlOutput += '<ul>';

    // تعریف ترتیب و ترجمه فارسی برای نمایش
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
        { key: 'matchType', label: 'نوع تطابق' },
        { key: 'confidence', label: 'درجه اطمینان' }
    ];

    displayOrder.forEach(item => {
        if (data[item.key]) { // فقط اگر کلید در داده وجود داشت، نمایش بده
            htmlOutput += `<li><strong>${item.label}:</strong> ${data[item.key]}</li>`;
        }
    });

    htmlOutput += '</ul>';
    resultsContainer.innerHTML = htmlOutput;
}

// تابع برای نمایش پیام خطا
function showError(message) {
    errorContainer.innerHTML = `<p>${message}</p>`;
    errorContainer.style.display = 'block';
    resultsContainer.innerHTML = ''; // پاک کردن نتایج قبلی اگر خطایی رخ داد
}
