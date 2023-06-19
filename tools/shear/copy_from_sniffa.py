import subprocess

def copy_files(system, specifics):
    SHEAR = ["3", "15", "60", "300", "1200"]
    REMOTE_PATH = f"salberti@sniffa:/scratch/salberti/MPMDP/production/shear_flow/atomistic/{system}/{specifics}/DATA"

    for shear_val in SHEAR:
        mom_file = f"mom_{shear_val}.txt"
        vel_file = f"vel_prof_{shear_val}_av.txt"
        remote_mom_file = f"{REMOTE_PATH}/{mom_file}"
        remote_vel_file = f"{REMOTE_PATH}/{vel_file}"           
        #subprocess.run(["scp", "-i", "/home/simon/.ssh/id_rsa", remote_mom_file, remote_vel_file, "."])
        #subprocess.run(["scp", "-i", "/c/Users/alber/.ssh/id_rsa", "-o", "IdentityFile=/c/Users/alber/.ssh/id_rsa", "-o", "LogLevel=DEBUG", remote_mom_file, remote_vel_file, "."])
        subprocess.run(["scp", "-i", "/c/Users/alber/.ssh/id_rsa", remote_mom_file, remote_vel_file, "."])

# Example usage:
system = "LJ2"
specifics = ""
#system = "hPF_force_massaging"
#specifics = "2x2x6"
copy_files(system, specifics)
