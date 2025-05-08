# Taxonomy App - Species Classification Search

A simple and engaging web application to search and display taxonomy (classification) information and images of living organisms using scientific APIs.

**✨ [View Live Demo](https://mynameiskaveh.github.io/taxonomy-app/) ✨** 
<!-- !!! Make sure this is your correct GitHub Pages link !!! -->

## 🌟 About The Project

This project was developed as a learning exercise and to showcase the ability to work with external APIs, build a serverless Python backend, and develop a frontend using HTML, CSS, and JavaScript. Users can enter the scientific or common name of an organism to view its classification hierarchy from kingdom to species, along with an image of the organism (if available).

## 🚀 Technologies Used

*   **Backend:**
    *   Python 3
    *   Flask (for creating a simple API)
    *   Serverless Functions hosted on Vercel
*   **Frontend:**
    *   HTML5
    *   CSS3
    *   JavaScript (Vanilla JS)
    *   Hosted on GitHub Pages
*   **External APIs:**
    *   [GBIF API (Global Biodiversity Information Facility)](https://www.gbif.org/developer/summary) - For retrieving taxonomic classification data.
    *   [Wikipedia API (via Python's `wikipedia` library)](https://pypi.org/project/wikipedia/) - For fetching images of organisms.
*   **Tools:**
    *   Git & GitHub (for version control and code hosting)

## ✨ Features

*   Search organisms by scientific name (and to some extent, common names supported by the APIs).
*   Display full taxonomic hierarchy (Kingdom, Phylum, Class, Order, Family, Genus, Species).
*   Show additional information like taxonomic rank, name status, and GBIF match confidence.
*   Display an image of the organism fetched from Wikipedia (if found).
*   Simple, clean, and user-friendly interface.
*   Responsive design for decent viewing on various devices.
*   Error handling and user-friendly feedback messages.

## 🛠️ Challenges & Learnings (Optional, but good to include)

*   Selecting reliable and suitable APIs for scientific data.
*   Handling different JSON response structures from various APIs.
*   Implementing logic to select the most relevant image from sources like Wikipedia, which may contain multiple or unrelated images.
*   Deploying a Python (Flask) application dificuldades serverlessly on Vercel.
*   Managing Cross-Origin Resource Sharing (CORS) requests between the frontend (GitHub Pages) and backend (Vercel).
*   Practicing coding and editing directly on the GitHub platform via its web interface.

## 💡 Future Enhancements (Optional)

*   Implement more advanced search functionalities (e.g., with more filters).
*   Allow users to choose from multiple results if a search query is ambiguous.
*   Integrate additional APIs for more comprehensive data or better image sourcing.
*   Further improve UI/UX aspects.

## 📝 License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
<!-- If you added a LICENSE.md file, this link will work. Otherwise, remove this line or just state "MIT License". -->

---

Developed by Kaveh
<!-- You can link to your GitHub profile here -->
