import streamlit as st
import pandas as pd
import re
from io import BytesIO

# ============================================================
# 1. Regex for rsID extraction
# ============================================================
RS_PATTERN = re.compile(r"\b(rs\d+)\b", re.IGNORECASE)

def extract_rs_id(text):
    if pd.isna(text):
        return None
    m = RS_PATTERN.search(str(text))
    return m.group(1).lower() if m else None


# ============================================================
# 2. Detect genotype column
# ============================================================
def find_genotype_column(df, user_id):
    user_id = user_id.strip()
    if not user_id:
        return None

    expected_1 = f"{user_id}.Plus/Minus Alleles"
    expected_2 = f"{user_id}-R.Plus/Minus Alleles"

    cols_lower = {c.lower(): c for c in df.columns}

    if expected_1.lower() in cols_lower:
        return cols_lower[expected_1.lower()]
    if expected_2.lower() in cols_lower:
        return cols_lower[expected_2.lower()]

    # Fallback: any column containing user ID
    for c in df.columns:
        if user_id.lower() in c.lower():
            return c

    return None


# ============================================================
# 3. Read variant list (small, OK to load fully)
# ============================================================
def read_variant_list(file_obj=None, text_input=""):
    if file_obj is not None:
        content = file_obj.getvalue().decode("utf-8", errors="replace")
    elif text_input.strip():
        content = text_input.strip()
    else:
        raise ValueError("No variant list provided.")

    tokens = re.split(r"[,\n\s]+", content)
    tokens = [x.strip() for x in tokens if x.strip()]

    rsids = set(
        filter(None, (extract_rs_id(t) for t in tokens))
    )

    if not rsids:
        raise ValueError("No recognizable rsIDs found (e.g., rs12345).")

    return rsids


# ============================================================
# 4. Process VCF-like file (CHUNKED READING)
# ============================================================
def process_vcf_file(uploaded_file, user_id, variant_set, chunksize=50000):

    filtered_chunks = []
    first_chunk = True

    for chunk in pd.read_csv(
        uploaded_file,
        sep="\t",
        dtype=str,
        chunksize=chunksize,
        low_memory=False
    ):

        # Detect genotype column only on first chunk
        if first_chunk:
            genotype_col = find_genotype_column(chunk, user_id)
            if genotype_col is None:
                raise ValueError(
                    f"Unable to find genotype column for user '{user_id}'.\n\n"
                    f"Available columns: {list(chunk.columns)}"
                )
            first_chunk = False

        # Rename required columns
        rename_map = {
            "Name": "Variant",
            "Chr": "Chr",
            "Position": "Pos",
            genotype_col: "Genotype"
        }
        chunk = chunk.rename(columns=rename_map)

        # Ensure required columns exist
        required = ["Variant", "Chr", "Pos", "Genotype"]
        if not all(k in chunk.columns for k in required):
            continue

        # Extract RSID
        chunk["NormalizedRSID"] = (
            chunk["Variant"]
            .str.extract(r"(?i)\b(rs\d+)\b", expand=False)
            .str.lower()
        )

        # Filter
        match_chunk = chunk[chunk["NormalizedRSID"].isin(variant_set)]

        # Keep only the four required columns
        match_chunk = match_chunk[["Variant", "Chr", "Pos", "Genotype"]]

        if len(match_chunk) > 0:
            filtered_chunks.append(match_chunk)

    # Combine all matching chunks
    if filtered_chunks:
        return pd.concat(filtered_chunks, ignore_index=True)
    else:
        return pd.DataFrame(columns=["Variant", "Chr", "Pos", "Genotype"])


# ============================================================
# 5. STREAMLIT UI
# ============================================================
st.set_page_config(page_title="Variant Data Extractor", layout="wide")
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

# MAIN UI INPUTS
uploaded_vcf = st.file_uploader(
    "1. Upload VCF-like File Tab-separated files",
    type=["tsv", "txt", "csv"]
)

user_id = st.text_input(
    "2. Enter User ID",
    value="IG1234"
).strip()

colA, colB = st.columns(2)

with colA:
    uploaded_variant_file = st.file_uploader(
        "Upload Variant List File (txt/csv)",
        type=["txt", "csv"]
    )

with colB:
    variant_list_raw = st.text_area(
        "Or paste rsID list here",
        height=100
    )


# ============================================================
# RUN
# ============================================================
if st.button("Extract", type="primary"):

    if uploaded_vcf is None:
        st.error("Please upload the VCF-like file.")
        st.stop()

    if not user_id:
        st.error("Please enter a valid user ID.")
        st.stop()

    try:
        with st.spinner("Reading Variant List..."):
            variant_set = read_variant_list(uploaded_variant_file, variant_list_raw)

        with st.spinner("Processing large VCF file (chunked)..."):
            df_filtered = process_vcf_file(uploaded_vcf, user_id, variant_set)

        if df_filtered.empty:
            st.warning("No matching variants found.")
        else:
            st.success(f"Found {len(df_filtered)} matching variants!")
            st.dataframe(df_filtered, use_container_width=True)

            tsv_bytes = df_filtered.to_csv(
                sep="\t", index=False
            ).encode("utf-8")

            st.download_button(
                label="📥 Download Extracted Variants",
                data=tsv_bytes,
                file_name=f"extracted_variants_{user_id}.tsv",
                mime="text/tab-separated-values"
            )

    except Exception as e:
        st.error(f"❌ Error: {e}")
