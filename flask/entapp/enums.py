"""
Detailed description.
"""








class TaskManagerState(enum.IntEnum):
    """
    Detailed description.
    """
    Idle = enum.auto()
    Error = enum.auto()
    RunningDownload = enum.auto()
    RunningConfig = enum.auto()
    RunningJob = enum.auto()








class TaskResult(enum.IntEnum):
    """
    Detailed description.
    """
    Finished = enum.auto()
    Error = enum.auto()
