"""

.. module:: taskqueue
    :synopsis:  multiprocessor stuffs

..  moduleauthor:: Ken Sugino <ken.sugino@gmail.com>

"""

import multiprocessing
from multiprocessing import TimeoutError
try:
    from Queue import Empty, Full
except:
    from queue import Empty, Full
import time
import traceback

import logging
logging.basicConfig(level=logging.DEBUG)
LOG = logging.getLogger(__name__)

# c.f. http://stackoverflow.com/questions/19924104/python-multiprocessing-handling-child-errors-in-parent

class Worker(multiprocessing.Process):
    
    def __init__(self, index, winfo, task_queue, result_queue):
        multiprocessing.Process.__init__(self)
        self.index = index
        self.winfo = winfo
        self.task_queue = task_queue
        self.result_queue = result_queue
        self.winfo[self.index] = {'status':'ready'}
        self.durs = []
        
    def set_info(self, **kw):
        d = self.winfo[self.index]
        d.update(kw)
        self.winfo[self.index] = d # need to trigger __set_item__
        
    def run(self):
        proc_name = self.name
        while True:
            wi = self.winfo[self.index]
            if wi.get('stop',False):
                print('{0}: Exiting (through stop)'.format(proc_name))
                self.set_info(status='exit')
                break
            try:
                next_task = self.task_queue.get(timeout=10)                
                if next_task is None:
                    # Poison pill means shutdown
                    maxdur = N.max(self.durs)
                    mindur = N.mean(self.durs)
                    numrun = len(self.durs)
                    print('{0}: Exiting (through shutdown) maxdur({1:.2f}) avgdur({2:.2f}) run({3})'.format(proc_name, maxdur, avgdur,numrun))
                    self.set_info(status='exit')
                    self.task_queue.task_done()
                    break
                print('{0}: starting {1}'.format(proc_name, next_task.name))
                stime = time.time()
                self.set_info(_stime=stime, status='running')
                try:
                    answer = next_task()
                    etime = time.time()
                    self.result_queue.put((next_task.name, answer))
                    elapsed = etime - stime
                    print('{0}: finished {1} ({2:.3} sec)'.format(proc_name, next_task.name, elapsed))
                    self.durs.append(elapsed)
                    self.set_info(_etime=etime, status='waiting', _stime=etime, elapsed=elapsed,
                        maxdur=N.max(self.durs), mindur=N.min(self.durs), avgdur=N.mean(self.durs))
                except Exception as e:
                    tb = traceback.format_exc()
                    print('{0} ({1}): error'.format(proc_name, next_task.name))
                    print(e)
                    print(tb)
                    self.set_info(error=str(e), trackback=str(tb),
                                 _etime=time.time(), status='waiting')
                finally:
                    self.task_queue.task_done()
            except Empty:
                pass
        return
    

class Task(object):
    def __init__(self, name, func, args=[], kwargs={}):
        self.func = func
        self.args = args
        self.kwargs = kwargs
        self.name = name

    def __call__(self):
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
        self.status = 'ready' # ready=>started=>stopped
        
    def start(self):
        if self.status=='ready':
            LOG.info('{1}: Starting {0} workers'.format(self.np, self.name))
            for w in self.workers:
                w.start()
            self.status = 'started'
        else:
            LOG.warninig('{0} already started'.format(self.name))
            
    def add_task(self, task):
        print('#Task({0}) added'.format(task.name))
        self.task_queue.put(task)
        
    def get_result(self, block=True, timeout=None):
        return self.result_queue.get(block,timeout)

    def set_info(self, i, **kw):
        d = self.winfo[i]
        d.update(kw)
        self.winfo[i] = d # need to trigger __set_item__

    def get_info(self, i):
        return self.winfo[i]

    def stop(self):
        if self.status=='started':
            for i in range(len(self.workers)):
                self.set_info(i, stop=True)
            self.status = 'stopped'
            for i in range(len(self.workers)):
                self.workers[i].join()
            LOG.info('{0} stopped'.format(self.name))
        else:
            LOG.warning('{0} not running (status: {1})'.format(self.name, self.status))
    
    def shutdown(self):
        if self.status == 'started':
            # send signals to workers
            for i in range(len(self.workers)):
                self.task_queue.put(None)
            # Wait for all of the tasks to finish
            self.task_queue.join()
            self.status = 'shutdown'
            LOG.info('{0} shutdown'.format(self.name))
        else:
            LOG.info('{0} shutdown for status ({1})'.format(self.name, self.status))

    def check_error(self, maxtime=None):
        for i in range(len(self.workers)):
            wi = self.get_info(i)
            err = wi.get('error', None)
            wname = self.workers[i].name
            if (err is not None):
                tb = wi['trackback']
                # stop whole thing
                print('STOPPING SERVER EXCEPTION in  {0}'.format(wname))
                print(err)
                print(tb)
                self.stop()
                return False
            if wi.get('_etime',0)<wi.get('_stime',0): # running
                elapsed = time.time()-wi['_stime']
            else:
                elapsed = 0
            if maxtime is not None:
                if (elapsed>maxtime):
                    print('STOPPING SERVER: elapsed {0}sec in  {1}'.format(elapsed, wname))
                    self.stop()
                    return False
        return True

    def __enter__(self):
        self.start()
        
    def __exit__(self, exc_type, exc_value, traceback):
        self.shutdown()


######        
def func1(a,b,c):
    time.sleep(0.5)
    return a+b+c

def func2(x):
    time.sleep(1)
    return x*x

def func3(x):
    time.sleep(0.1)
    raise RuntimeError('test error message')

def loop():
    num_workers = 2 
    server = Server(num_workers)
    # Enqueue initial jobs
    num_jobs = 5
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
        running = True
        while server.check_error(4):
            print('####### loop head #######')
            try:
                result = server.get_result(timeout=5) # block until result come in
            except Empty:
                result = None
                print('queue empty server: {0}'.format(server.name))
            if result is not None:
                print('Result: name:{0}, answer:{1}'.format(*result))
                wi = server.winfo
                # for k,v in wi.items():
                #     print('  {0:.2f}: elapsed'.format(time.time()-v['_stime']))
                if result[0].startswith('phase1_'):
                    tname = 'phase2_{0}'.format(result[1])
                    args = (result[1],)
                    t = Task(tname, func2, args)
                    server.add_task(t)
                    p1[result[0]] = result[1]
                    # if result[1] == 3:
                    #     print('!!!!!!!!!!!!!!!!!!!!!!!  putting in func3')
                    #     server.add_task(Task('ERRORTASK',func3,(1,)))                    
                if result[0].startswith('phase2_'):
                    p2[result[0]] = result[1]
                    if result[1] == 9:
                        print('!!!!!!!!!!!!!!!!!!!!!!!  putting in func3')
                        server.add_task(Task('ERRORTASK',func3,(1,)))                    
                # if len(p2)==num_jobs:
                #     break
        print('Exit Loop')
    print('Done')        
    print(p1)
    print(p2)