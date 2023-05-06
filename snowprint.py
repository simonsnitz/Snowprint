from src.Create_Alignment import create_alignment
from src.Create_Regulators import create_regulators
from src.Create_Operons import create_operons
from src.Update_Associations import update_associations
from src.Create_Operators import create_operators
from src.Analyze_Result import pull_operator, write_frontend_json

from alive_progress import alive_bar

import sys





    # Perform the Snowprint workflow
def predict_operator(acc, **kwargs):

    with alive_bar(5, dual_line=True, title='Snowprint') as bar:
      
        bar.text = f'-> Collecting homologs for {acc}'
        create_alignment(acc)
        bar()

        bar.text = f'-> Collecting metadata for {acc} homologs'
        create_regulators(acc)
        bar()

        bar.text = f'-> Collecting genetic context for {acc} homologs'
        create_operons(acc)
        bar()

        bar.text = f'-> Extracting operators for {acc} homologs'
        update_associations(acc)
        bar()

        if 'known_operator' in kwargs:
            bar.text = f'-> Aligning homologous operator and generating metrics for {acc}'
            create_operators(acc, known_operator=kwargs.get('known_operator'))
            bar()
        else:
            bar.text = f'-> Aligning operators and generating metrics for {acc}'
            create_operators(acc)
            bar()






    # Perform the Snowprint workflow AND update the data.json file for frontend display
def snowprint(acc, **kwargs):

    if 'known_operator' in kwargs:
        predict_operator(acc, known_operator=kwargs.get('known_operator'))
        data = pull_operator(acc)
        write_frontend_json(data)
    else:
        predict_operator(acc)
        data = pull_operator(acc)
        write_frontend_json(data)






if __name__ == "__main__":

    acc = sys.argv[1]
    snowprint(acc)
