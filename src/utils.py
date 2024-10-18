def pretty_print_predictions(predictions: list[dict[str, float]]) -> str:
    """
    Pretty print the ligand class predictions in a formatted way.

    :param predictions: A list of dictionaries containing 'Class' as the ligand class name
                    and 'Probability' as the associated probability.
    :return: A formatted HTML string containing the ligand class predictions
    """
    # Start building the HTML string
    html_output = """
        <b>Ligand Class Predictions</b><br>
        <b>(Click on ligand group name to see the full list of ligands)</b><br><br>
        <table border="1" cellpadding="5" cellspacing="0">
            <tr>
                <th><b>Ligand Class</b></th>
                <th><b>Probability</b></th>
            </tr>
        """

    # Add rows for each prediction, with the ligand class as a clickable link
    base_url = "https://checkmyblob.bioreproducibility.org/server/ligands/#"
    for prediction in predictions:
        ligand_class = prediction['Class']
        probability = prediction['Probability']
        link = f'<a href="{base_url}{ligand_class}" target="_blank">{ligand_class}</a>'
        html_output += f"<tr><td>{link}</td><td>{probability:.3f}</td></tr>"

    # Close the HTML table
    html_output += "</table>"

    return html_output
