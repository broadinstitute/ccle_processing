import tempfile

import dalmatian as dm
from genepy.utils.helper import createFoldersFor
from depmapomics.config import CACHE_PATH

class OmicsPipeline():
    def __init__(self, workspace, tmp_path=None, verbose=True):
        assert CACHE_PATH[-1] == '/'
        createFoldersFor(CACHE_PATH)
        self.verbose = verbose
        self.workspace = workspace
        self.sample_delete_protection = True
        if tmp_path is None:
            self.tmp_path = tempfile.mkdtemp(dir='/tmp/')
        else:
            self.tmp_path = tmp_path
        self.printv('fetching the workspace')
        self.workspace_manager = dm.WorkspaceManager(workspace)
        self.printv('fetching sample sets from the workspace')
        self.samplesets = self.workspace_manager.get_sample_sets()
        self.printv('fetching samples from the workspace')
        self.samples = self.workspace_manager.get_samples()
        self.printv('fetching participants from the workspace')
        self.participants = self.workspace_manager.get_participants()

    def printv(self, text):
        if self.verbose:
            print(text)

    def delete_samples_on_terra(self):
        if self.sample_delete_protection:
            print('Doing nothing. self.sample_delete_protection needs to be set to False manually for this to work.')
        else:
            self.sample_delete_protection = True
            self.printv('deleting samples info on terra')
            self.workspace_manager.disable_hound()
            self.workspace_manager.delete_sample(self.samples.index)
            self.workspace_manager.delete_sample_set(self.sample_sets.index)
            self.workspace_manager.delete_participant(self.participants.index)
