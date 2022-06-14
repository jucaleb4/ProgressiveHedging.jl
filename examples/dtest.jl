using Pkg
Pkg.activate("..")
using Distributed 

for i in Distributed.workers()
    id, pid, host = fetch(@spawnat i (myid(), getpid(), gethostname()))
    println("my id: ",id,", ", "hostname: ", host)
end
