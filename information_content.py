def calculate_ic_depth_breadth(ic_file):
    """
    Calculates IC depth, IC breadth, and normalized IC breadth.

    Parameters:
    - ic_file (str): Path to the file with IC values per protein (tab-separated, one protein per line).

    Returns:
    - dict: A dictionary with keys 'ic_depth', 'ic_breadth', 'normalized_ic_breadth'.
    """
    total_ic = 0.0
    total_annotations = 0
    num_proteins = 0

    with open(ic_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            ic_values = [float(x) for x in line.split('\t') if x.strip()]
            if ic_values:
                total_ic += sum(ic_values)
                total_annotations += len(ic_values)
                num_proteins += 1

    ic_depth = total_ic / total_annotations if total_annotations else 0.0
    ic_breadth = total_ic
    normalized_ic_breadth = total_ic / num_proteins if num_proteins else 0.0

    return {
        'ic_depth': ic_depth,
        'ic_breadth': ic_breadth,
        'normalized_ic_breadth': normalized_ic_breadth
    }

