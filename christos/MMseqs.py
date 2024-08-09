import subprocess
import os
import sys
from pathlib import Path
import shutil

class MMseqs2API:

    def __init__(self, mmseqs2_path='/home/simon.herman/.local/bin/mmseqs/bin/mmseqs', threads=1, cleanup=False):

        Path(f"{os.getcwd()}/mmseqs_temp").mkdir(parents=True, exist_ok=True)

        self.mmseqs2_path = mmseqs2_path
        self.threads = threads
        self.cleanup = cleanup
        self.dir = f"{os.getcwd()}/mmseqs_temp"

    def _run_command(self, command, log_file):

        print(f"Running command: {command}")

        with open(log_file, 'w') as log:

            process = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)

            for stdout_line in iter(process.stdout.readline, ""):
                log.write(stdout_line)

            for stderr_line in iter(process.stderr.readline, ""):
                print(stderr_line, end="", file=sys.stderr)
                log.write(stderr_line)


            process.stdout.close()
            process.stderr.close()
            return_code = process.wait()

            if return_code != 0:
                raise Exception(f"Command '{command}' failed with return code {return_code}")

    def createdb(self, fasta, log_file='createDB.log'):

        path = Path(fasta)
        if not path.exists():
            raise FileNotFoundError(f"FASTA file {fasta} does not exist.")
        
        command = f"{self.mmseqs2_path} createdb {path} {self.dir}/db"
        self._run_command(command, log_file)

    def cluster(self, coverage, identity, cov_mode, log_file='clustering.log'):

        command = (
            f"{self.mmseqs2_path} cluster {self.dir}/db {self.dir}/db_clu {self.dir}/tmp -c {coverage} --min-seq-id {identity} --cov-mode {cov_mode} --cluster-mode 2 --cluster-reassign -v 3 --threads {self.threads}"
        )
        self._run_command(command, log_file)

    def createtsv(self, tsv_file='result_seq_clu.tsv', log_file='createTsv.log'):
        command = f"{self.mmseqs2_path} createtsv {self.dir}/db {self.dir}/db {self.dir}/db_clu {self.dir}/{tsv_file}"
        self._run_command(command, log_file)

        if self.cleanup:
            self._cleanup_intermediate_files(tsv_file)
    
    def _cleanup_intermediate_files(self, tsv_file):
        files_to_keep = {f"{self.dir}/{tsv_file}"}
        for file in Path(self.dir).iterdir():
            if str(file) not in files_to_keep:
                try:
                    if file.is_file():
                        file.unlink()
                    elif file.is_dir():
                        shutil.rmtree(file)
                except Exception as e:
                    print(f"Error deleting file {file}: {e}")


# Example usage
if __name__ == '__main__':


    mmseqs2_api = MMseqs2API(threads=40, cleanup=True)
    
    # Example file lists
    fasta_file_single = '/store/EQUIPES/BIM/MEMBERS/paul.roginski/OLD/ORFPRED/new_data/disordered/DisProt_release_2023_06_CUT.faa'
    
    # Running the commands with a single file
    mmseqs2_api.createdb(fasta_file_single)
    mmseqs2_api.cluster(coverage=0.7, identity=0.4, cov_mode=0)
    mmseqs2_api.createtsv()
