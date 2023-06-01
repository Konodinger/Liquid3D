import docker
client = docker.from_env()

# run the container
command = "ls -l"
print(client.containers.run("barthpaleologue/openvdb_bridge", command))