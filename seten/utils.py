from os import path, mkdir


def filename(filepath):
    """
    Obtain filename from the filepath
    """
    return path.splitext(path.basename(filepath))[0]

def output_results(output, results):
    """
    Outputs gene set enrichment results with their gene counts and p-values
    """
    # output directory exists?
    if not path.exists(path.dirname(output)):
        mkdir(path.dirname(output))

    # any result available
    if results:
        with open(output, 'w') as f:
            header = results[0].keys()
            writer = csv.DictWriter(f, delimiter='\t', fieldnames=header)
            writer.writeheader()
            for result in results:
                if result:
                    writer.writerow(result)
