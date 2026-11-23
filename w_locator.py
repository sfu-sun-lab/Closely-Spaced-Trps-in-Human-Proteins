"""
Code follow the flow:
    Compute W position and count
    Compute intra-W distance (Output "position_w_output.xlsx")
    Compute W frequencies (Output "frequency_w_output.csv")
    Compute r gap standard deviation and Z-score for the first 10 gaps (Output "final_w_output.xlsx")

Usage:
    for example, run from the command line:
    python w_locator.py
    or, for Mac/Linux:
    python3 w_locator.py

Dependencies: pandas, numpy, openpyxl (pip install pandas numpy openpyxl)
"""
import ast
from pathlib import Path
from openpyxl.styles import Font, Alignment

import numpy as np
import pandas as pd
import re

from xlsxwriter.worksheet import Worksheet

SEQUENCE_INPUT = "sequenceInput_2024_05_14.xlsx"
TM_ORIENTATION_INPUT = "W_and_TM_input.csv"
POSITION_W_OUTPUT = "position_w_output.xlsx"
FREQ_OUTPUT_CSV = "frequency_w_output.csv"
FINAL_OUTPUT = "final_w_output.xlsx"

ENTRY_COLUMN = "Entry"
POSITION_COLUMN = "Sequence W position"
COUNT_COLUMN = "Sequence W count"
LENGTH_COLUMN = "Length"
PROTEIN_NAMES_COLUMN = "Protein names"
GENE_NAMES_COLUMN = "Gene Names"

FINAL_PROTEIN_NAMES_COLUMN = "Protein name"
FINAL_GENE_NAMES_COLUMN = "Gene symbol"
FINAL_COUNT_COLUMN = "Trp count"

_EMPTY = np.array([], dtype=np.int64)

### COMPUTE W FREQUENCY CODE ###

def parse_positions(raw):
    """Parse one position-list cell into sorted, de-duplicated positions.

    Returns (positions, malformed). Missing values, blank strings, invalid
    literals, and non-integer elements are treated as malformed rows with
    no valid positions. Duplicates are removed before distances are
    calculated so they can never create a spurious zero-distance.
    """
    if raw is None or (isinstance(raw, float) and np.isnan(raw)):
        return _EMPTY, True
    text = raw.strip() if isinstance(raw, str) else raw
    if text == "":
        return _EMPTY, True
    try:
        parsed = ast.literal_eval(text) if isinstance(text, str) else text
    except (ValueError, SyntaxError, TypeError):
        return _EMPTY, True

    if isinstance(parsed, int):
        parsed = [parsed]
    elif not isinstance(parsed, (list, tuple)):
        return _EMPTY, True

    try:
        ints = [int(x) for x in parsed]
    except (TypeError, ValueError):
        return _EMPTY, True

    return np.unique(np.asarray(ints, dtype=np.int64)), False


def calculate_frequencies(positions):
    """Per-entry (distance value, count) pairs via np.diff + np.unique."""
    if positions.size < 2:
        return _EMPTY, _EMPTY
    distances = np.diff(positions)
    return np.unique(distances, return_counts=True)


def build_output(entries, row_vals, row_counts):
    """Assemble the Entry + frequency-columns table in one vectorized pass.

    Columns are created only for distances observed anywhere in the
    dataset, and the whole matrix is filled with a single fancy-index
    assignment rather than per-row/per-column insertion.
    """
    lengths = np.fromiter((v.size for v in row_vals), dtype=np.int64, count=len(entries))
    if lengths.sum() == 0:
        return pd.DataFrame({ENTRY_COLUMN: entries})

    all_vals = np.concatenate(row_vals)
    all_counts = np.concatenate(row_counts)
    distances = np.unique(all_vals)  # sorted, observed-only

    row_index = np.repeat(np.arange(len(entries)), lengths)
    col_index = np.searchsorted(distances, all_vals)

    max_count = int(all_counts.max())
    dtype = next(
        d for d in (np.uint8, np.uint16, np.uint32, np.uint64)
        if max_count <= np.iinfo(d).max
    )
    matrix = np.zeros((len(entries), distances.size), dtype=dtype)
    matrix[row_index, col_index] = all_counts

    columns = [f"Inter-Trp freq {int(d)}" for d in distances]
    output = pd.DataFrame(matrix, columns=columns)
    output.insert(0, ENTRY_COLUMN, entries)
    return output

def compute_w_frequency(df: pd.DataFrame) -> pd.DataFrame:
    print("Begin Computing W Frequency Table")
    for column in (ENTRY_COLUMN, POSITION_COLUMN):
        if column not in df.columns:
            raise ValueError(f"Required column not found in data: {column!r}")

    entries = df[ENTRY_COLUMN].tolist()

    # Information list to keep
    lengths = df[LENGTH_COLUMN].to_numpy()
    w_counts = df[COUNT_COLUMN].to_numpy()
    proteins = df[PROTEIN_NAMES_COLUMN].to_numpy()
    genes = df[GENE_NAMES_COLUMN].to_numpy()

    row_vals, row_counts = [], []
    malformed_rows = 0
    total_gaps = 0

    for raw in df[POSITION_COLUMN]:
        positions, malformed = parse_positions(raw)
        malformed_rows += malformed
        vals, counts = calculate_frequencies(positions)
        total_gaps += int(counts.sum())
        row_vals.append(vals)
        row_counts.append(counts)

    output = build_output(entries, row_vals, row_counts)
    # Drop empty frequency column
    freq_cols = [c for c in output.columns if c.startswith("Inter-Trp freq ")]
    empty_freqs = output[freq_cols].ne(0).any(axis=1).to_numpy()
    output = output[empty_freqs]

    output.to_csv(FREQ_OUTPUT_CSV, index=False)

    max_distance = max((int(c.rsplit(" ", 1)[1]) for c in freq_cols), default=0)

    print("Processing complete.")
    print(f"Input rows: {len(df)}")
    print(f"Total inter-Trp gaps: {total_gaps}")
    print(f"Unique observed distances: {len(freq_cols)}")
    print(f"Maximum inter-Trp distance: {max_distance}")
    print(f"Malformed rows: {malformed_rows}")
    print(f"Output shape: {output.shape[0]} rows x {output.shape[1]} columns")
    print(f"Saved to: {FREQ_OUTPUT_CSV}")

    output.insert(1, FINAL_PROTEIN_NAMES_COLUMN, proteins[empty_freqs])
    output.insert(2, FINAL_GENE_NAMES_COLUMN, genes[empty_freqs])
    output.insert(3, LENGTH_COLUMN, lengths[empty_freqs])
    output.insert(4, FINAL_COUNT_COLUMN, w_counts[empty_freqs])

    return output

### COMPUTE W GAP STD AND Z-SCORE CODE ###

def process_w_frequency(data: pd.DataFrame):
    print("Begin Processing W Frequency Table")
    # Keep first ten gap
    data_columns = [ENTRY_COLUMN, FINAL_PROTEIN_NAMES_COLUMN, FINAL_GENE_NAMES_COLUMN, LENGTH_COLUMN, FINAL_COUNT_COLUMN]
    distances = [d for d in range(1, 11) if f"Inter-Trp freq {d}" in data.columns]
    columns_to_compute = [f"Inter-Trp freq {d}" for d in distances]

    # Drop row if Inter-Trp freq up to 10 is empty
    data = data[data[columns_to_compute].ne(0).any(axis=1)]

    data_columns.extend(columns_to_compute)
    output = data[data_columns]
    N = output["Trp count"].to_numpy()
    L = output[LENGTH_COLUMN].to_numpy()
    prob = N / L
    result_dict = {}
    for distance in distances:

        x = output[f"Inter-Trp freq {distance}"].to_numpy()
        mean = (N - 1) * prob * np.power((1 - prob), distance - 1)
        q = prob * np.power((1 - prob), distance - 1)
        std = np.sqrt((N - 1) * q * (1 - q))
        z = np.divide((x - mean), std, out=np.zeros_like(std), where=std != 0)
        result_dict[distance] = {"std": std, "z": z}

    for distance in distances:
        output.insert(len(output.columns), f"Inter-Trp standard deviation {distance}", result_dict[distance]['std'])
    for distance in distances:
        output.insert(len(output.columns), f"Inter-Trp Z-score {distance}", result_dict[distance]['z'])

    # Sort data by W gap-1 frequency
    output = output.sort_values([f"Inter-Trp freq {distances[0]}", ENTRY_COLUMN], ascending=False)
    output.set_index(ENTRY_COLUMN, inplace=True)

    # Create multi-level columns
    final_columns = [(FINAL_PROTEIN_NAMES_COLUMN, ''), (FINAL_GENE_NAMES_COLUMN, ''), (LENGTH_COLUMN, ''), (FINAL_COUNT_COLUMN, '')]
    for distance in distances:
        final_columns.append(("Frequency", f"W-gap {distance}"))
    for distance in distances:
        final_columns.append(("Standard deviation", f"W-gap {distance}"))
    for distance in distances:
        final_columns.append(("Z score", f"W-gap {distance}"))
    output.columns = pd.MultiIndex.from_tuples(final_columns, names=[ENTRY_COLUMN, ''])

    with pd.ExcelWriter(FINAL_OUTPUT, engine='openpyxl') as writer:
        sheet_name = 'Sheet1'
        num_levels = output.columns.nlevels
        output.to_excel(writer, sheet_name=sheet_name)

        worksheet = writer.sheets[sheet_name]

        # Delete the annoying empty row when use multi-level columns
        worksheet.delete_rows(num_levels+1)

        # Merge single level header
        for col in range(1, worksheet.max_column + 1):
            top_header = worksheet.cell(row=1, column=col).value
            sub_headers = [worksheet.cell(row=row, column=col).value for row in range(2, num_levels+1)]
            if top_header is not None and all(str(value).strip() == '' or value is None for value in sub_headers):
                worksheet.merge_cells(start_row=1, start_column=col, end_row=num_levels, end_column=col)

        # Update Header alignment
        for row in range(1, num_levels + 1):
            for col in range(1, worksheet.max_column + 1):
                    worksheet.cell(row, col).alignment = Alignment(horizontal='center', vertical='center')

        #Remove bold font for entry
        for row in range(num_levels+1, len(output.index)+1):
            worksheet.cell(row, 1).font = Font(name='Calibri')
            worksheet.cell(row, 1).alignment = Alignment(horizontal='center', vertical='bottom')

    print(f"Saved to: {FINAL_OUTPUT}")

### COMPUTE W POSITION CODE ###

def compute_w_position() -> pd.DataFrame:
    print("Begin Computing W Position")
    def get_integer_value(newdata):
        string_list = re.findall('\d+', newdata)
        numerical_list = list(map(int, string_list))
        return numerical_list

    input_path = Path(SEQUENCE_INPUT)
    if not input_path.exists():
        raise FileNotFoundError(f"Required input file not found: {SEQUENCE_INPUT}")

    data : pd.DataFrame = pd.read_excel(input_path)

    # Calculating length of sequences
    position = pd.Series([
        ''.join(f',{j + 1}' for j, char in enumerate(seq) if char == 'W')
        for seq in data['Sequence']
    ])

    data[POSITION_COLUMN] = position.apply(get_integer_value)
    count = pd.Series(data[POSITION_COLUMN])
    length = count.str.len()
    data[COUNT_COLUMN] = length

    return data

### EXTRACT INTRA-W CODE ###

class IntraTMR:
    def __init__(self, position: int, tm_left: int, tm_right: int):
        self.position = position
        self.tm_left = tm_left
        self.tm_right = tm_right

    def get_left_distance(self, is_absolute = True) -> int:
        distance = self.position - self.tm_left
        if is_absolute:
            distance = abs(distance)
        return distance

    def get_right_distance(self) -> int:
        return self.tm_right - self.position

    def get_closet_distance(self) -> int:
        return min(self.get_left_distance(), self.get_right_distance())

    def __str__(self):
        return f"W{self.position}({self.get_closet_distance()})"

def extract_tm_position(series):
    tm_position = []
    tm_literal = []
    if series[' Num of TMRs'] != 0:
        topologies = series[' Orientation'].strip()
        tm_pattern = re.compile(r'(\d+)-(\d+)')
        split_topologies = re.split('[io]', topologies)
        for tm_extract in split_topologies:
            m = tm_pattern.match(tm_extract)
            if not m:
                continue
            tm_position.append(f'{m[1]}-{m[2]}')
            tm_literal.append([int(m[1]), int(m[2])])
    series['TM positions'] = tm_position
    series['TM literal'] = tm_literal
    return series

def locate_intra_w(series):
    tm_positions = series['TM literal']
    w_positions = series[POSITION_COLUMN]
    intra_distance_info = {}
    intra_tmr = []
    intra_tmr_position = []
    if tm_positions and w_positions:
        for tm in tm_positions:
            intra_distance_info[f'{tm}'] = {'W': [], 'distance': []}
            for w in w_positions:
                if w in intra_tmr_position:
                    continue
                if tm[0] <= w <= tm[1]:
                    intra_tmr_position.append(w)
                    intra_tmr.append(IntraTMR(w, tm[0], tm[1]))
    series['Intra-W'] = ", ".join(str(x) for x in intra_tmr)

    return series

def extract_tryptophan_distance(pos_data):
    print("Begin Computing Intra-W Distance")
    orientation_data = pd.read_csv(TM_ORIENTATION_INPUT)
    orientation_data.rename(columns={'Name': ENTRY_COLUMN}, inplace=True)
    orientation_data = orientation_data[[ENTRY_COLUMN, " Orientation", " Num of TMRs"]]

    df = pos_data.merge(orientation_data, on=ENTRY_COLUMN)

    df = df.apply(extract_tm_position, result_type='expand', axis=1)
    df = df.apply(locate_intra_w, result_type='expand', axis=1)

    # Drop empty column
    df = df.loc[df['Intra-W'] != '', :]

    # Drop Intermediate columns
    output_df = df.drop([" Orientation", "TM literal"], axis=1)
    output_df.to_excel(POSITION_W_OUTPUT, index=False, engine='openpyxl')
    print(f"Saved to: {POSITION_W_OUTPUT}")

def main():
    pos_data = compute_w_position()
    extract_tryptophan_distance(pos_data)
    #pos_data = pd.read_csv("Total_W.csv")
    freq_data = compute_w_frequency(pos_data)
    process_w_frequency(freq_data)

if __name__ == "__main__":
    main()