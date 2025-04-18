import argparse
import logging

logging.basicConfig(level=logging.INFO, format="[%(levelname)s]: %(message)s")


def parse_arguments() -> argparse.Namespace:
    """Parse command line arguments for the script.

    Returns:
        argparse.Namespace: Parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description="Generate a vcoord file for the supercomputer Fugaku.")

    parser.add_argument("-f", "--file_name", type=str, default="./vcoord_file", help="Name of the output file (default: ./vcoord_file)")
    parser.add_argument("-x", "--x", type=int, default=8, help="Number of X dimensions (default: 8)")
    parser.add_argument("-y", "--y", type=int, default=8, help="Number of Y dimensions (default: 8)")
    parser.add_argument("-z", "--z", type=int, default=8, help="Number of Z dimensions (default: 8)")
    parser.add_argument("-p", "--proc", type=int, default=4, help="Number of processes (default: 4)")
    parser.add_argument("-t", "--thread", type=int, default=12, help="Number of threads per process (default: 12)")

    return parser.parse_args()


def generate_vcoord_file(
        file_name: str,
        x: int,
        y: int,
        z: int,
        proc: int,
        thread: int,
    ) -> None:
    """Generate a vcoord file of the supercomputer Fugaku.

    This configuration aims to minimize the communication distance between adjacent computational regions.

    Args:
        file_name (str): Name of the output file.
        x (int): Number of X dimensions.
        y (int): Number of Y dimensions.
        z (int): Number of Z dimensions.
        proc (int): Number of processes.
        thread (int): Number of threads per process.
    """
    # Check number of cores to use
    core = proc * thread
    assert core > 0, f"Invalid core count: {proc}x{thread}={core}. The number of cores must be positive."
    assert core <= 48, f"Core count exceeds Fugaku's limit: {proc}x{thread}={core}. Maximum is 48 cores."
    if core == 48:
        logging.info("Configuration uses all 48 available cores.")
    elif core < 48:
        logging.warning(f"Configuration does not use all available cores. Current setting uses {core} cores.")

    # Show setting info
    print()
    logging.info(f"Node configuration: {x}x{y}x{z}")
    logging.info(f"Number of processes: {proc}")
    logging.info(f"Number of threads: {thread}")

    # Generate base lines
    lines = []
    for k in range(z):
        for j in range(y):
            for i in range(x):
                line = "(" + str(i) + "," + str(j) + "," + str(k) + ")"
                line += " core=" + str(thread)
                lines.append(line)

    # Reverse Y direction
    for idx in range(0, len(lines), x*y):
        if idx % (2*x*y) == 0:
            continue
        else:
            lines[idx:idx+x*y] = lines[idx:idx+x*y][::-1]
            for jdx in range (idx, idx+x*y, x):
                lines[jdx:jdx+x] = lines[jdx:jdx+x][::-1]

    # Multiple lines
    lines = [s for s in lines for _ in range(proc)]

    # Write file
    print()
    try:
        with open(file_name, "w") as f:
            f.write("\n".join(lines))
        logging.info(f"Output file: {file_name}")
        logging.info(f"Job script description: \"node={x}x{y}x{z}:torus:strict\"")
    except IOError as e:
        logging.error(f"Failed to write to file {file_name}: {e}")


def main():
    args = parse_arguments()
    generate_vcoord_file(args.file_name, args.x, args.y, args.z, args.proc, args.thread)


if __name__ == "__main__":
    main()
