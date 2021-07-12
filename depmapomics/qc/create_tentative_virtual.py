from taigapy import TaigaClient
from depmapomics.qc.config import TENTATIVE_VIRTUAL, TAIGA_IDS_LATEST


def get_taiga_id_with_version(taiga_id):
    tc = TaigaClient()
    latest_version = max([int(x['name']) for x in tc.get_dataset_metadata(taiga_id)['versions']])
    taiga_id_with_version = '{}.{}'.format(taiga_id, latest_version)
    return taiga_id_with_version


def get_taiga_ids_list(taiga_id_dict):
    taiga_ids_list = []
    for taiga_id, pairs in taiga_id_dict.items():
        taiga_id_with_version = get_taiga_id_with_version(taiga_id)
        for pair in pairs:
            taiga_ids_list.append({
                'name': pair[0],
                'taiga_id': '{}/{}'.format(taiga_id_with_version, pair[1])
            })
    return taiga_ids_list


def update_tentative_virtual():
    tc = TaigaClient()

    # Uncomment if the tentative virtual dataset doesn't exist
    #  with open('/tmp/null.csv', 'w') as f:
    #     pass

    # new_dataset_id = tc.create_dataset(
    #     "tentative virtual",
    #     dataset_description="this is a temporary virtual dataset for omics QCs", # optional (but recommended)
    #     upload_files=[
    #         {
    #             'path': '/tmp/null.csv',
    #             'format': 'Raw'
    #         }
    #     ],
    #     folder_id="29b48277847443edb7d84b9b457de124", # optional, will default to your home folder if not provided
    # )


    taiga_ids_list = get_taiga_ids_list(TAIGA_IDS_LATEST)

    new_dataset_id = tc.update_dataset(
        TENTATIVE_VIRTUAL['name'],
        changes_description="add all the virtual data",
        add_taiga_ids=taiga_ids_list,
        add_all_existing_files = False
    )


if __name__ == "__main__":
    update_tentative_virtual()
