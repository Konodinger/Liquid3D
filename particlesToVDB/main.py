import docker
import os
client = docker.from_env()

def wrapCommand(command):
    return "bash -c '" + command + "'"

DOCKER_IMAGE_NAME = "barthpaleologue/openvdb_bridge"

# Set the container directory path where you want to mount the host directory
CONTAINER_DIRECTORY = "/root/Liquid3D/openvdbCppBridge/mountedHostDirectory"

HOST_DIRECTORY = os.getcwd()

# Create a dictionary to configure the volume mount
volume_config = {
    HOST_DIRECTORY: {
        'bind': CONTAINER_DIRECTORY,
        'mode': 'rw'
    }
}


def execute_commands(commands: list[str]):
    container = client.containers.run(image=DOCKER_IMAGE_NAME, command='/bin/sh', tty=True, detach=True,
                                    # mounts current directory to /root/Liquid3D/openvdbCppBridge/build/ in container
                                    volumes=volume_config)

    for COMMAND in commands:
        exit_code, result = container.exec_run(wrapCommand(COMMAND), stream=True, privileged=True, stdout=True, stderr=True)
        if exit_code:
            print("Error running command: " + COMMAND)
            print("Exit code: ", exit_code)
            break

        for line in result:
            print(line.decode('utf-8'))

    print("Stopping container")
    container.stop()
    print("Container stopped")
    print("Removing container")
    container.remove()
    print("Container removed")


def convert_to_vdb(file_name):
    # check if file exists
    if not os.path.isfile(file_name):
        print("File does not exist:", file_name)
        return

    file_path_in_container = os.path.join(CONTAINER_DIRECTORY, file_name)
    execute_commands([
        f"/root/Liquid3D/openvdbCppBridge/build/OpenVdbBridge -f {file_path_in_container}",
        f'ls /results',
        f"cp -r /results {CONTAINER_DIRECTORY}/",
    ])

convert_to_vdb("./liquidPointCloud.txt")