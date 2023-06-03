import docker
import os
client = docker.from_env()

DOCKER_IMAGE_NAME = "barthpaleologue/openvdb_bridge"

# Set the container directory path where you want to mount the host directory
CONTAINER_DIRECTORY = "/root/Liquid3D/openvdbCppBridge/build/"

HOST_DIRECTORY = os.getcwd()

# Create a dictionary to configure the volume mount
volume_config = {
    HOST_DIRECTORY: {
        'bind': CONTAINER_DIRECTORY,
        'mode': 'rw'
    }
}

COMMAND = "/root/Liquid3D/openvdbCppBridge/build/OpenVdbBridge -f /root/Liquid3D/openvdbCppBridge/liquidPointCloud_002.txt"

COMMAND2 = "ls -a /root/Liquid3D/openvdbCppBridge/build/"

COMMAND3 = "echo 'hello world' > /root/Liquid3D/openvdbCppBridge/build/test.txt"

COMMAND4 = "cat /root/Liquid3D/openvdbCppBridge/build/test.txt"

"""lines = client.containers.run(DOCKER_IMAGE_NAME, COMMAND, volumes={os.getcwd(): {'bind': '/tmp/', 'mode': 'rw'}}, stream=True)                                                                         

for line in lines:                                                                                  
    print(line)"""


container = client.containers.run(image=DOCKER_IMAGE_NAME, command='/bin/sh', tty=True, detach=True, 
                                  # mounts current directory to /root/Liquid3D/openvdbCppBridge/build/ in container
                                  volumes=volume_config)

result = container.exec_run(COMMAND)
print(result.output.decode('utf-8'))
result = container.exec_run(COMMAND2)
print(result.output.decode('utf-8'))
result = container.exec_run(COMMAND3)
print(result.output.decode('utf-8'))
result = container.exec_run(COMMAND4)
print(result.output.decode('utf-8'))

print("Stopping container")
container.stop()
print("Container stopped")
print("Removing container")
container.remove()
print("Container removed")