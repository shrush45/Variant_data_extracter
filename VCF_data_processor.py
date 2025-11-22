import streamlit as st
import pandas as pd
import re
from io import StringIO, BytesIO

# -------------------------------------------------------------
# 1. Extract rsID (vectorized & clean)
# -------------------------------------------------------------
def extract_rs_id(text):
    if pd.isna(text):
        return None
    match = re.search(r"\b(rs\d+)\b", str(text), flags=re.IGNORECASE)
    return match.group(1).lower() if match else None

# -------------------------------------------------------------
# 2. Genotype Column Detection (with contains user_id feature)
# -------------------------------------------------------------
def find_genotype_column(df, user_id):
    user_id = user_id.strip()
    if not user_id:
        return None

    expected_1 = f"{user_id}.Plus/Minus Alleles"
    expected_2 = f"{user_id}-R.Plus/Minus Alleles"

    cols_lower = {c.lower(): c for c in df.columns}

    # Try exact matches first
    if expected_1.lower() in cols_lower:
        return cols_lower[expected_1.lower()]
    if expected_2.lower() in cols_lower:
        return cols_lower[expected_2.lower()]

    # Fallback — any column containing user id
    for c in df.columns:
        if user_id.lower() in c.lower():
            return c

    return None

# -------------------------------------------------------------
# 3. Process VCF File
# -------------------------------------------------------------
def process_vcf_file(uploaded_file, user_id):
    try:
        raw_bytes = uploaded_file.getvalue()
        raw_text = raw_bytes.decode("utf-8", errors='replace')
    except Exception as e:
        raise ValueError(f"Error reading file: {e}")

    df = pd.read_csv(StringIO(raw_text), sep="\t", dtype=str)

    # Detect genotype column
    genotype_col = find_genotype_column(df, user_id)
    if genotype_col is None:
        raise ValueError(
            f"Unable to find genotype column for user '{user_id}'.\n\n"
            "Expected patterns:\n"
            f"  - {user_id}.Plus/Minus Alleles\n"
            f"  - {user_id}-R.Plus/Minus Alleles\n"
            "Or *any* column containing the user ID.\n"
            f"Available columns: {list(df.columns)}"
        )

    rename_map = {
        "Name": "Variant",
        "Chr": "Chr",
        "Position": "Pos",
        genotype_col: "Genotype"
    }

    df = df.rename(columns=rename_map)

    required = ["Variant", "Chr", "Pos", "Genotype"]
    missing = [c for c in required if c not in df.columns]

    if missing:
        raise ValueError(f"Missing columns after processing: {missing}")

    return df[required]


# -------------------------------------------------------------
# 4. Read Variant List
# -------------------------------------------------------------
def read_variant_list(file_obj=None, text_input=""):
    if file_obj is not None:
        try:
            content = file_obj.getvalue().decode("utf-8", errors="replace")
        except:
            raise ValueError("Unable to read uploaded variant list file.")
    elif text_input.strip():
        content = text_input.strip()
    else:
        raise ValueError("No variant list provided.")

    tokens = re.split(r"[,\n\s]+", content)
    tokens = [x.strip() for x in tokens if x.strip()]

    normalized = set(filter(None, (extract_rs_id(t) for t in tokens)))

    if not normalized:
        raise ValueError("No recognizable rsIDs found (e.g., rs12345).")

    return normalized


# -------------------------------------------------------------
# 5. Filter Variants
# -------------------------------------------------------------
def filter_variants(df, variant_set):
    df["NormalizedRSID"] = (
        df["Variant"]
        .str.extract(r"(?i)\b(rs\d+)\b", expand=False)
        .str.lower()
    )
    out = df[df["NormalizedRSID"].isin(variant_set)].copy()
    out = out.drop(columns=["NormalizedRSID"], errors="ignore")
    return out


# -------------------------------------------------------------
# 6. STREAMLIT UI SETUP
# -------------------------------------------------------------
st.set_page_config(
    page_title="Variant Data Extractor",
    layout="wide",
    initial_sidebar_state="expanded"
)

st.title("🧬 Variant Data Extractor")


# -------------------------------------------------------------------
# SIDEBAR (same style as your original script)
# -------------------------------------------------------------------
st.sidebar.header("Instructions")
st.sidebar.markdown("""
1. **Upload VCF-like file** (tab-separated).
2. **Enter User ID** (must match genotype column pattern).
3. **Provide rsID list** by file upload or text box.
4. Click **Run** to extract matching variants.
5. Use **Reset** to clear all inputs.
""")

st.sidebar.info(
    "Script automatically detects genotype columns using "
    "`USERID`, `USERID-R`, or any column containing the user ID."
)

st.sidebar.markdown("---")


# -------------------------------------------------------------------
# MAIN BODY UI
# -------------------------------------------------------------------
uploaded_vcf = st.file_uploader(
    "1. Upload VCF-like Tab-Separated File",
    type=["txt", "tsv", "csv"],
    key="vcf_file"
)

user_id = st.text_input("2. Enter User ID", value="IG1234", key="user_id_input").strip()

st.subheader("Provide Variant List")

colA, colB = st.columns(2)

with colA:
    uploaded_variant_file = st.file_uploader(
        "Upload Variants File",
        type=["txt", "csv"],
        key="variant_file"
    )

with colB:
    variant_list_raw = st.text_area(
        "OR paste rsID list here",
        value="rs12976533\nrs267333\nrs34670133",
        height=100,
        key="variant_textarea"
    )

st.markdown("---")


# -------------------------------------------------------------------
# RUN + RESET BUTTONS (side-by-side)
# -------------------------------------------------------------------
col_run, col_reset = st.columns(2)

run_clicked = col_run.button("🚀 Run", type="primary")
reset_clicked = col_reset.button("🔄 Reset", type="primary")

if reset_clicked:
    # Option 1 (recommended) — clear everything stored in session_state, then rerun
    st.session_state.clear()
    st.rerun()
# -------------------------------------------------------------------
# PROCESSING PIPELINE
# -------------------------------------------------------------------
if run_clicked:

    if uploaded_vcf is None:
        st.error("Please upload the main VCF-like file.")
        st.stop()

    if not user_id:
        st.error("Please enter a valid User ID.")
        st.stop()

    try:
        with st.spinner("Processing VCF file..."):
            df_processed = process_vcf_file(uploaded_vcf, user_id)

        with st.spinner("Reading variant list..."):
            variant_set = read_variant_list(uploaded_variant_file, variant_list_raw)

        with st.spinner("Filtering variants..."):
            df_filtered = filter_variants(df_processed, variant_set)

        if df_filtered.empty:
            st.warning("No matching variants found.")
        else:
            st.success(f"Found **{len(df_filtered)}** matching variants!")
            st.dataframe(df_filtered, use_container_width=True)

            tsv_bytes = df_filtered.to_csv(
                sep="\t", index=False
            ).encode("utf-8")

            st.download_button(
                label="📥 Download Extracted data",
                data=tsv_bytes,
                file_name=f"extracted_variants_{user_id}.tsv",
                mime="text/tab-separated-values"
            )

    except Exception as e:
        st.error(f"❌ Error: {e}")
