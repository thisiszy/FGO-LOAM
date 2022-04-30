import os
import sys
import subprocess
import time
    
def detect_active():
    val = os.popen("rosnode list")
    for line in val.readlines():
        if("/rosbag_play" in line.split()[0]):
            return True
    return False
def get_ape(seq, filename):
    # val = subprocess.check_output("evo_ape kitti tmp.txt "+filename+" -a")
    val = os.popen("evo_ape kitti /data/dataset/kitti/GT/"+seq+".txt "+filename+" -a")
    print("------------------------")
    for line in val.readlines():
        tmp = line.split()
        print(tmp)
        if(len(tmp) == 2 and tmp[0] == "rmse"):
            return tmp[1]

if __name__ == '__main__':
    seq = sys.argv[1]
    output_folder = sys.argv[2]
    dataset = sys.argv[3]
    output_folder += seq + '/'
    dataset += seq
    print(dataset)
    print(output_folder)

    len_dic = {"00":4541, "01":1101, "02":4661, "03":801, "04":271,
               "05":2761, "06":1101, "07":1101, "08":4071, "09":1591, "10":1201}

    fail = 0
    start_time = time.time()
    while(True):
        result_file_name = output_folder + "dt_kitti.txt"
        sleep_est = len_dic[seq] * 30 / 1000
        print("seq: ",seq," start")
        launch_proc = subprocess.Popen(["roslaunch", "lcloam", "lcloam_eval.launch", "output_folder:="+output_folder, "bag_path:="+dataset+".bag"], shell=False)
        time.sleep(10)
        while(detect_active()):
            print("rosbag has played for {0} seconds...".format(str(time.time()-start_time)))
            time.sleep(3)
        time.sleep(sleep_est + fail * 100)
        print("----seq end-----")
        launch_proc.terminate()
        time.sleep(10)
        val = os.popen("wc -l "+result_file_name)
        tmp = val.readline()
        if(tmp.split()[0] == str(len_dic[seq])):
            exit()
        else:
            fail += 1
            print("\033[1;32mERROR once\033[0m")

        if(fail > 1):
            print("\033[0;31;40mERROR 2 times\033[0m")
            exit()
