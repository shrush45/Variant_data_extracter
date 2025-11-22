#  Variant Data Extractor (Streamlit App)

This is a **Streamlit** web application designed to efficiently process VCF-like (tab-separated) genotype files, automatically identify the user's genotype column, and filter the data to extract specific variants based on a list of **rsIDs**.

##  Features

* **Robust Genotype Column Detection:** Automatically finds the genotype column using patterns like `USERID.Plus/Minus Alleles`, `USERID-R.Plus/Minus Alleles`, or any column containing the User ID.
* **Vectorized rsID Extraction:** Uses efficient string operations to normalize and extract rsIDs (e.g., `rs12345`) from both the VCF-like file and the user-provided variant list.
* **Flexible Variant Input:** Accepts the list of target rsIDs via a file upload or a direct text area input.
* **Interactive UI:** Provides a clean, user-friendly interface built with Streamlit for easy file upload, parameter setting, and result display.
* **Downloadable Output:** Allows users to download the filtered results as a tab-separated file (`.tsv`).

##  Prerequisites

To run this application, you need **Python** and the following libraries:

* `streamlit`
* `pandas`

##  Installation

1.  **Clone the repository** (or download `VCF_data_processor.py`):
    ```bash
    git clone https://github.com/shrush45/Variant_data_extracter
    cd Variant_data_extracter
    ```

2.  **Install the required libraries:**
    ```bash
    pip install streamlit pandas
    ```

##  How to Run

Execute the script using the Streamlit CLI:

```bash
streamlit run VCF_data_processor.py