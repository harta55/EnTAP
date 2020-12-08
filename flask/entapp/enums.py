"""
Detailed description.
"""
import enum








class TaskManagerState(enum.IntEnum):
    """
    Detailed description.
    """
    Idle = enum.auto()
    Error = enum.auto()
    Running = enum.auto()








class TaskResult(enum.IntEnum):
    """
    Detailed description.
    """
    Running = enum.auto()
    Finished = enum.auto()
    Error = enum.auto()
