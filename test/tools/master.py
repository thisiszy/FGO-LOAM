import sys
import subprocess
import signal
if __name__ == '__main__':
    output_folder = sys.argv[1]
    dataset_folder = sys.argv[2]
    print(output_folder)
    print(dataset_folder)

    len_dic = {"00":4541, "01":1101, "02":4661, "03":801, "04":271,
                "05":2761, "06":1101, "07":1101, "08":4071, "09":1591, "10":1201}
    try:
        for seq in ["00", "01", "02", "03", "04", "05", "06", "07", "08", "09", "10"]:
            #  #  #  #  #  #  #  #  #  #  #  # 
            process = subprocess.Popen(["python3", "run_record_subproc.py", seq, output_folder, dataset_folder, "> rubish.txt"])
            process.wait()
            #  #  #  #  #  #  #  #  #  #  #  # 
        print("All test done!")
    except KeyboardInterrupt:
        process.send_signal(signal.SIGINT)
        exit()