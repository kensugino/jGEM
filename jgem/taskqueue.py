import multiprocessing
import time

class Worker(multiprocessing.Process):
    
    def __init__(self, index, winfo, task_queue, result_queue):
        multiprocessing.Process.__init__(self)
        self.index = index
        self.winfo = winfo
        self.task_queue = task_queue
        self.result_queue = result_queue
        self.winfo[self.index] = {}
        
    def set_info(self, **kw):
        d = self.winfo[self.index]
        d.update(kw)
        self.winfo[self.index] = d # need to trigger __set_item__
        
    def run(self):
        proc_name = self.name
        while True:
            next_task = self.task_queue.get()
            if next_task is None:
                # Poison pill means shutdown
                print('{0}: Exiting'.format(proc_name))
                self.set_info(status='exit')
                self.task_queue.task_done()
                break
            print('{0}: starting {1}'.format(proc_name, next_task.name))
            self.set_info(_stime=time.time(), status='running')
            answer = next_task()
            self.task_queue.task_done()
            self.result_queue.put((next_task.name, answer))
            print('{0}: finished {1}'.format(proc_name, next_task.name))
            self.set_info(_etime=time.time(), status='waiting')
        return
    

class Task(object):
    def __init__(self, name, func, args=[], kwargs={}):
        self.func = func
        self.args = args
        self.kwargs = kwargs
        self.name = name

    def __call__(self):
        #print('Task {0} start'.format(self.name))
        return self.func(*self.args)
        
class Server(object):
    
    def __init__(self, np=2, name='Server'):
        self.np = np
        self.name = name
        self.result_queue = rq = multiprocessing.Queue()
        self.task_queue = tq = multiprocessing.JoinableQueue()
        self.manager = multiprocessing.Manager()
        self.winfo = wi = self.manager.dict()
        self.workers = [ Worker(i, wi, tq, rq) for i in range(np)]
        
    def start(self):
        print('{1}: Starting {0} workers'.format(self.np, self.name))
        for w in self.workers:
            w.start()
            
    def add_task(self, task):
        self.task_queue.put(task)
        
    def get_result(self, block=True, timeout=None):
        return self.result_queue.get(block,timeout)
    
    def shutdown(self):
        for i in range(len(self.workers)):
            self.task_queue.put(None)
        # Wait for all of the tasks to finish
        self.task_queue.join()
        print('{0} shutdown'.format(self.name))
        
    def __enter__(self):
        self.start()
        
    def __exit__(self, exc_type, exc_value, traceback):
        self.shutdown()
#################### 

class Worker2(multiprocessing.Process):
    
    def __init__(self, task_queue):
        multiprocessing.Process.__init__(self)
        self.task_queue = task_queue

    def run(self):
        proc_name = self.name
        while True:
            next_task = self.task_queue.get()
            if next_task is None:
                # Poison pill means shutdown
                print('{0}: Exiting'.format(proc_name))
                self.task_queue.task_done()
                break
            print('{0}: {1}'.format(proc_name, next_task.name))
            next_task()
            self.task_queue.task_done()
        return
    

class Task2(object):
    
    def __init__(self, name, result_queue, func, args=[], kwargs={}):
        self.func = func
        self.args = args
        self.name = name
        self.result_queue = result_queue
        self.kwargs = kwargs
        
    def __call__(self):
        print('Task {0} start'.format(self.name))
        rslt = self.func(*self.args, **self.kwargs)
        self.result_queue.put((self.name, rslt))
        
#
class Server2(object):
    
    def __init__(self, np=2, name='Server'):
        self.np = np
        self.name = name
        self.task_queue = tq = multiprocessing.JoinableQueue()
        self.workers = [ Worker2(tq) for i in range(np)]
        
    def start(self):
        print('Server {1}: Starting {0} workers'.format(self.np, self.name))
        for w in self.workers:
            w.start()
            
    def add_task(self, task):
        self.task_queue.put(task)
            
    def shutdown(self):
        for i in range(len(self.workers)):
            self.task_queue.put(None)
        # Wait for all of the tasks to finish
        self.task_queue.join()
        print('Server {0} shutdown'.format(self.name))
        
    def __enter__(self):
        self.start()
        
    def __exit__(self, exc_type, exc_value, traceback):
        self.shutdown()
        

######      
def loop2():
    server = Server(2)
    manager = multiprocessing.Manager()
    rq1 = manager.Queue()
    #rq2 = multiprocessing.Queue()
    # Wait for results
    p1 = {}
    p2 = {}
    with server:
        # Enqueue initial jobs
        num_jobs = 10
        for i in range(num_jobs):
            tname = 'phase1_{0}'.format(i)
            t = Task(tname, rq1, func1, (i,i,i))
            server.add_task(t)
        while True:
            name,ans = rq1.get()
            print('Result: name:{0}, answer:{1}'.format(name,ans))
            if name.startswith('phase1_'):
                tname = 'phase2_{0}'.format(ans)
                t = Task(tname, rq1, func2, kwargs={'x':ans} )
                server.add_task(t)
                p1[name] = ans
            if name.startswith('phase2_'):
                p2[name] = ans
            if len(p2)==num_jobs:
                break
        print('Exit Loop')
    print('Done')


######        
def func1(a,b,c):
    time.sleep(0.2)
    return a+b+c

def func2(x):
    time.sleep(0.1)
    return x*x

def loop():
    server = Server(2)
    # Enqueue initial jobs
    num_jobs = 10
    for i in range(num_jobs):
        tname = 'phase1_{0}'.format(i)
        args = (i,i,i)
        t = Task(tname, func1, args)
        server.add_task(t)
    # Wait for results
    print(server.winfo)
    p1 = {}
    p2 = {}
    with server:
        while True:
            result = server.get_result()
            print('Result: name:{0}, answer:{1}'.format(*result))
            wi = server.winfo
            for k,v in wi.items():
                print('  {0:.2f}: elapsed'.format(time.time()-v['_stime']))
            if result[0].startswith('phase1_'):
                tname = 'phase2_{0}'.format(result[1])
                args = (result[1],)
                t = Task(tname, func2, args)
                server.add_task(t)
                p1[result[0]] = result[1]
            if result[0].startswith('phase2_'):
                p2[result[0]] = result[1]
            if len(p2)==num_jobs:
                break
        print('Exit Loop')
    print('Done')        
    print(p1)
    print(p2)