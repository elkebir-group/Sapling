import sys

def handle_args():
    boolean = {"m"}
    parameter = {"a":0.9,"t":5,"l":-1,"m":False}
    waiting = None
    for arg in sys.argv[1:]:
        if waiting is None:
            if arg[0]!="-":
                print("unrecongized options",file=sys.stderr)
                exit(0)
            waiting = arg[1:]
            if waiting in boolean:
                parameter[waiting] = True
                waiting = None
        else :
            parameter[waiting] = arg
            waiting = None
    return parameter

def handle_arg_expand():
    boolean = {"m"}
    parameter = {"a":0.9,"m":False}
    waiting = None
    for arg in sys.argv[1:]:
        if waiting is None:
            if arg[0]!="-":
                print("unrecongized options",file=sys.stderr)
                exit(0)
            waiting = arg[1:]
            if waiting in boolean:
                parameter[waiting] = True
                waiting = None
        else :
            parameter[waiting] = arg
            waiting = None
    return parameter