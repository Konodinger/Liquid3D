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


COMMANDS = [
    "/root/Liquid3D/openvdbCppBridge/build/OpenVdbBridge -f /root/Liquid3D/openvdbCppBridge/liquidPointCloud_002.txt",
    f"mkdir {CONTAINER_DIRECTORY}/docker_results",
    #f"cp /root/Liquid3D/openvdbCppBridge/build/OpenVdbBridge {CONTAINER_DIRECTORY}/docker_results/",
    #f"cp /root/Liquid3D/openvdbCppBridge/liquidPointCloud_002.txt {CONTAINER_DIRECTORY}/docker_results/",
    #f"{CONTAINER_DIRECTORY}/docker_results/OpenVdbBridge -f {CONTAINER_DIRECTORY}/docker_results/liquidPointCloud_002.txt",
    f'ls /results',
    f"cp -r /results {CONTAINER_DIRECTORY}/docker_results/",
    f"ls /root/Liquid3D/openvdbCppBridge/build/",
    f"touch {CONTAINER_DIRECTORY}/docker_results/proof.txt",
    f"echo 'proof' >> {CONTAINER_DIRECTORY}/docker_results/proof.txt"
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