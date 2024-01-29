from skaha.session import Session

# if True record memory and cpu usage as a function of time using psrecord and scalene
memcpu_logging = True

# job will execute this command with argument galaxyID
if memcpu_logging:
    cmd = "/arc/projects/mauve/products/scripts/log_make_products.sh"
else:
    cmd = "/arc/projects/mauve/products/scripts/make_products.sh"

# select container image
image = "images.canfar.net/skaha/base:latest"

# galaxies to iterate over
galaxy_list = ["IC3392"]  # , "NGC4383", "NGC4064", "NGC4694"]

# loop over galaxies and launch jobs
for galaxy in galaxy_list:
    galaxy = galaxy.upper().replace(" ", "")

    args = galaxy
    print("Launching job")
    print("Command:\n{} {}".format(cmd, args))
    print("")

    # smaller cubes require less RAM
    cores = 16
    if galaxy in ["NGC4501"]:
        ram = 192  # GB

    elif galaxy not in ["NGC4694"]:
        ram = 64  # GB
    else:
        ram = 32  # GB

    # launch headless session on the science platform
    session = Session()
    session_id = session.create(
        name=galaxy,
        image=image,
        cores=cores,
        ram=ram,
        kind="headless",
        cmd=cmd,
        args=args,
        env={"sessiontype": "headless"},
    )

    print("Sesion ID: {}".format(session_id[0]))
    print("")
