"""Shared helpers used by more than one analysis stage.

Modules here hold only code that was verified to behave identically at every
call site. Stage scripts keep their own orchestration; nothing in this package
reads or writes a result table.
"""
