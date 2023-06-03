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


def convert_to_vdb():
    COMMANDS = [
        "/root/Liquid3D/openvdbCppBridge/build/OpenVdbBridge -f /root/Liquid3D/openvdbCppBridge/liquidPointCloud_002.txt",
        f'ls /results',
        f"cp -r /results {CONTAINER_DIRECTORY}/",
    ]

    container = client.containers.run(image=DOCKER_IMAGE_NAME, command='/bin/sh', tty=True, detach=True,
                                    # mounts current directory to /root/Liquid3D/openvdbCppBridge/build/ in container
                                    volumes=volume_config)

    for COMMAND in COMMANDS:
        result = container.exec_run(wrapCommand(COMMAND), stream=True, privileged=True)
        for line in result.output:
            print(line.decode('utf-8'))

    print("Stopping container")
    container.stop()
    print("Container stopped")
    print("Removing container")
    container.remove()
    print("Container removed")

convert_to_vdb()