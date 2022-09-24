import os
import pandas as pd

col_names: list[str] = ['tool', 'tool_setup']


def get_performance_metrics_per_row(row: pd.Series):
    return pd.concat(
        [
            pd.read_csv(row['file'], sep='\t').reset_index(drop=True),
            pd.DataFrame(row).transpose().reset_index(drop=True)
        ],
        axis=1
    )


def get_performance_metrics(paths: pd.Series) -> pd.DataFrame:
    ret: list[pd.DataFrame] = []
    df = paths.str.split('/', expand=True).drop(columns=[0])
    dataset = df[df.shape[1]].str.extract(r'(.*)\.tsv$').rename(columns={0: 'dataset'})
    df.drop(columns=[df.shape[1]], inplace=True)
    df.columns = [col_names[i] for i in range(df.shape[1])]
    df = pd.concat(
        [df, dataset, paths],
        axis=1
    )
    for i in range(df.shape[0]):
        ret.append(
            get_performance_metrics_per_row(df.loc[i])
        )
    return pd.concat(ret, axis=0)


def main():
    ret: list[pd.DataFrame] = []
    for (root, dirs, files) in os.walk('.', topdown=True):
        df = pd.DataFrame(
            {
                'file_names': pd.Series(files),
                'is_target': pd.Series(files).str.match(r'^\d+\.tsv$')
            }
        ).assign(
            root=root
        )

        df['file'] = df.apply(lambda x: os.path.join(x['root'], x['file_names']), axis=1)

        if df['is_target'].any():
            ret.append(get_performance_metrics(df['file'][df['is_target']]))
    pd.concat(ret, axis=0).reset_index(drop=True).to_csv('file.tsv', sep='\t', index=False)


if __name__ == "__main__":
    main()
