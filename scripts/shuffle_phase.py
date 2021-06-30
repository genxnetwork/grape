import vcfpy
import random


if __name__ == '__main__':
    input_path = 'test_data/test.vcf.gz'
    output_path = 'test_data/shuffled.vcf.gz'
    reader = vcfpy.Reader.from_path(input_path)

    writer = vcfpy.Writer.from_path(output_path, header=reader.header)
    for record in reader:

        if not record.is_snv():
            continue

        new_calls = []
        for i, call in enumerate(record.calls):
            data = call.data.get('GT')
            if random.random() > 0.0:
                new_data = data[::-1]
                call.set_genotype(new_data)
                new_calls.append(call)

        record.update_calls(new_calls)
