from .cli import main

if __name__ == "__main__":
    # required on Windows: multiprocessing uses `spawn` and re-imports __main__
    import multiprocessing as mp

    mp.freeze_support()
    main()
