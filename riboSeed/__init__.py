from ._version import __version__
# this may be uneeded, but we get "cannot connect to X server " errors
# due to the import order.  we need to get mpl.use set early, and whats earlier
# than an init?
try:
    import matplotlib as mpl
    mpl.use('Agg')
except:
    pass
