import os
import glob
def find_none_hf(index_file = 'index.txt'):

    with open(index_file, "r") as f:
        lines = f.readlines()


    for line in lines[3:]:
        parts = line.strip().split("\t")
        if len(parts) >= 4:
            Hf = parts[-1]
            if Hf == "None":
                return True
    return False
def main():

    for mhfr_file in sorted(glob.glob("*.mhfr"), key=lambda x: int(x.split(".")[0])):
            if os.path.isdir(mhfr_file):
                try:
                    ind_path = os.path.join(mhfr_file, "index.txt")
                    missing = find_none_hf(ind_path)
                    if  missing:
                        print(f'Missing Hf in {mhfr_file}')
                except Exception as e:
                        print(f'Error in {mhfr_file}')
            else:
                print(f'Directory issue in {mhfr_file}')
def main_cli():
    main()        
        
if __name__ == "__main__":
    main()
    


